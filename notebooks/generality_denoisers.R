#CLR vs raw |Spearman| readout on a fixed imputed matrix, swept across denoisers: raw, kNN, SVD, ALRA, MAGIC
suppressPackageStartupMessages({
  library(Rmagic); library(Matrix); library(dplyr); library(tidyr); library(ggplot2)
  library(data.table); library(irlba); library(rsvd); library(RANN)
})
stopifnot(Rmagic::pymagic_is_available())
data_dir <- normalizePath(file.path("..","data"), mustWork=FALSE)
out_dir  <- normalizePath(file.path("..","notebooks"), mustWork=FALSE)
fig_dir  <- normalizePath(file.path("..","figures"), mustWork=FALSE)
NPCA<-20L; K<-30L; SEEDS<-101:105; RANK<-20L

#TF -> target lists
tf_mouse<-list(
  Gata1=c("Klf1","Hba-a1","Hba-a2","Hbb-b1","Hbb-b2","Aqp1","Slc4a1","Gypa","Epor"),
  Gata2=c("Gata1","Klf1","Itga2b"),Klf1=c("Hba-a1","Hba-a2","Hbb-b1","Hbb-b2","Aqp1","Slc4a1"),
  Sfpi1=c("Cebpa","Mpo","Csf1r","Csf2ra","Csf2rb"),Cebpa=c("Mpo","Elane","Prtn3","Ctsg","Csf3r"),
  Cebpe=c("Mpo","Elane","Prtn3","Ctsg"),Gfi1=c("Csf1r","Elane"),Irf8=c("Cxcr4"),Tal1=c("Gata1","Klf1","Gata2"))
tf_tusi<-tf_mouse; names(tf_tusi)[names(tf_tusi)=="Sfpi1"]<-"Spi1"
tf_tusi<-lapply(tf_tusi,function(v){v[v=="Hbb-b1"]<-"Hbb-bs";v[v=="Hbb-b2"]<-"Hbb-bt";v})
tf_human<-list(
  PAX5=c("CD19","MS4A1","CD79A","CD79B","BLNK","BANK1","BCL11A"),TBX21=c("IFNG","GZMB","CXCR3","GNLY","PRF1","NKG7"),
  EOMES=c("IFNG","GZMA","PRF1","GNLY","NKG7","KLRD1"),FOXP3=c("IL2RA","CTLA4","TNFRSF18"),
  SPI1=c("CD14","LYZ","CSF1R","FCGR3A","CST3","S100A8","S100A9"),IRF8=c("CCR7","FLT3","CST3","HLA-DRA"),
  CEBPB=c("CD14","LYZ","S100A8","S100A9"),LEF1=c("SELL","CCR7","TCF7"),KLF2=c("SELL","S1PR1","CCR7"),ID3=c("SELL","CCR7"))

#helpers
build_panel<-function(M,tf,min_lib,gene_pos,seed,n_cells=1500L,n_genes=2000L){
  ls<-rowSums(M);keep<-which(ls>=min_lib&ls<=quantile(ls,0.99));set.seed(seed)
  ci<-sample(keep,min(n_cells,length(keep)));cn<-M[ci,,drop=FALSE]
  mk<-intersect(unique(c(names(tf),unlist(tf))),colnames(cn))
  gp<-colMeans(cn>0);cn<-cn[,gp>=gene_pos|colnames(cn)%in%mk,drop=FALSE]
  gv<-apply(cn,2,var);hv<-order(gv,decreasing=TRUE)[seq_len(min(n_genes,length(gv)))]
  X<-cn[,union(colnames(cn)[hv],mk),drop=FALSE];storage.mode(X)<-"integer";X}
norm_sqrt<-function(X){s<-rowSums(X);s[s==0]<-1;Xn<-X/s*mean(s);sqrt(Xn)}
fast_auc<-function(scores,labels){r<-rank(scores);np<-sum(labels==1);nn<-sum(labels==0)
  if(np==0||nn==0)return(NA_real_);(sum(r[labels==1])-np*(np+1)/2)/(np*nn)}
make_pairs<-function(tf,gs,seed=7L){set.seed(seed);rows<-list()
  for(f in names(tf)){if(!(f%in%gs))next;tg<-intersect(tf[[f]],gs);if(!length(tg))next
    pool<-setdiff(gs,c(f,tg));ng<-sample(pool,min(length(pool),5L*length(tg)))
    rows[[f]]<-rbind(data.frame(tf=f,gene=tg,label=1L),data.frame(tf=f,gene=ng,label=0L))}
  do.call(rbind,rows)}
score_raw <- function(C,pr) abs(mapply(function(a,b)C[a,b],pr$tf,pr$gene))
score_clr <- function(C,pr){ mu<-rowMeans(C); sg<-apply(C,1,sd); sg[sg==0]<-1
  zi<-pmax(0,(mapply(function(a,b)C[a,b],pr$tf,pr$gene)-mu[pr$tf])/sg[pr$tf])
  zj<-pmax(0,(mapply(function(a,b)C[a,b],pr$tf,pr$gene)-mu[pr$gene])/sg[pr$gene])
  sqrt(zi^2+zj^2) }

#denoisers, each returns a cells x genes matrix with gene colnames
den_raw  <- function(X) { Y<-norm_sqrt(X); colnames(Y)<-colnames(X); Y }
den_knn  <- function(X, k=K, npca=NPCA){              #neighbour pooling in PCA space
  set.seed(1L); Xn<-norm_sqrt(X); pc<-irlba::prcomp_irlba(Xn,n=npca,center=TRUE,scale.=FALSE)$x
  idx<-RANN::nn2(pc,k=k)$nn.idx
  Y<-matrix(0,nrow(Xn),ncol(Xn)); for(i in seq_len(nrow(Xn))) Y[i,]<-colMeans(Xn[idx[i,],,drop=FALSE])
  colnames(Y)<-colnames(X); Y }
den_svd  <- function(X, k=RANK){                       #low-rank truncated SVD
  set.seed(1L); Xn<-norm_sqrt(X); sv<-rsvd::rsvd(Xn,k=k); Y<-sv$u%*%diag(sv$d,k,k)%*%t(sv$v)
  colnames(Y)<-colnames(X); Y }
den_alra <- function(X, k=RANK){                       #ALRA, Linderman et al 2022
  set.seed(1L); ls<-rowSums(X); ls[ls==0]<-1; A<-log1p(sweep(X,1,ls,"/")*median(ls))
  sv<-rsvd::rsvd(A,k=k); Ak<-sv$u%*%diag(sv$d,k,k)%*%t(sv$v)
  for(j in seq_len(ncol(Ak))){
    thr<-abs(min(Ak[,j])); col<-Ak[,j]; col[col<thr]<-0
    onz<-A[,j][A[,j]>0]; nz<-col>0
    if(sum(nz)>1 && length(onz)>1){ cm<-mean(col[nz]); cs<-sd(col[nz])
      if(!is.na(cs)&&cs>0) col[nz]<-(col[nz]-cm)/cs*sd(onz)+mean(onz); col[col<0]<-0 }
    Ak[,j]<-col }
  colnames(Ak)<-colnames(X); Ak }
den_magic<- function(X){ Y<-as.matrix(Rmagic::magic(X,t="auto",knn=K,npca=NPCA,seed=1L,verbose=0)$result)
  colnames(Y)<-colnames(X); Y }
DENOISERS<-list(raw=den_raw, knn=den_knn, svd=den_svd, alra=den_alra, magic=den_magic)

eval_denoisers<-function(X,tf){
  gs<-colnames(X); pr<-make_pairs(tf,gs); out<-list()
  for(dn in names(DENOISERS)){
    Y<-tryCatch(DENOISERS[[dn]](X), error=function(e){message("  ",dn," failed: ",conditionMessage(e)); NULL})
    if(is.null(Y)) next
    C<-suppressWarnings(cor(Y,method="spearman")); C[is.na(C)]<-0
    out[[dn]]<-data.frame(denoiser=dn,
      raw_auc=fast_auc(score_raw(C,pr),pr$label),
      clr_auc=fast_auc(score_clr(C,pr),pr$label)) }
  do.call(rbind,out) }

#loaders
load_paul<-function(){cf<-data.table::fread(file.path(data_dir,"GSE72857_umitab.txt"),sep="\t",header=TRUE,data.table=FALSE);gn<-cf[[1]];cf<-as.matrix(cf[,-1]);rownames(cf)<-gn;M<-t(cf);storage.mode(M)<-"integer";M}
load_pbmc<-function(){d<-file.path(data_dir,"filtered_gene_bc_matrices","hg19");m<-Matrix::readMM(file.path(d,"matrix.mtx"));gt<-read.table(file.path(d,"genes.tsv"),sep="\t");rownames(m)<-make.unique(gt$V2);m<-as.matrix(t(m));storage.mode(m)<-"integer";m}
load_tm<-function(){d<-file.path(data_dir,"droplet","Marrow-10X_P7_2");m<-Matrix::readMM(file.path(d,"matrix.mtx"));gt<-read.table(file.path(d,"genes.tsv"),sep="\t");rownames(m)<-make.unique(gt$V1);m<-as.matrix(t(m));storage.mode(m)<-"integer";m}
load_tusi<-function(){csv<-file.path(data_dir,"GSM2388072_basal_bone_marrow.raw_umifm_counts.csv");df<-data.table::fread(csv,sep=",",header=TRUE,data.table=FALSE);meta<-c("cell_id","barcode","library_id","seq_run_id","pass_filter");M<-as.matrix(df[,setdiff(colnames(df),meta)]);rownames(M)<-df$cell_id;storage.mode(M)<-"integer";M}
specs<-list(
  Paul=list(load=load_paul,tf=tf_mouse,min_lib=500L,gene_pos=0.05),
  PBMC=list(load=load_pbmc,tf=tf_human,min_lib=500L,gene_pos=0.03),
  TM  =list(load=load_tm,  tf=tf_mouse,min_lib=500L,gene_pos=0.03),
  Tusi=list(load=load_tusi,tf=tf_tusi, min_lib=1000L,gene_pos=0.05))

#optional env filters for a quick subset run
.only <- Sys.getenv("GEN_ONLY");  if(nzchar(.only))  specs<-specs[strsplit(.only,",")[[1]]]
.seeds<- Sys.getenv("GEN_SEEDS"); if(nzchar(.seeds)) SEEDS<-as.integer(strsplit(.seeds,",")[[1]])

#run
rows<-list()
for(nm in names(specs)){ sp<-specs[[nm]]; cat(sprintf("\n=== %s ===\n",nm)); M<-sp$load()
  for(s in SEEDS){ X<-build_panel(M,sp$tf,sp$min_lib,sp$gene_pos,seed=s)
    e<-eval_denoisers(X,sp$tf); e$dataset<-nm; e$seed<-s; rows[[paste(nm,s)]]<-e
    cat(sprintf("  seed %d | ",s))
    for(i in seq_len(nrow(e))) cat(sprintf("%s raw=%.3f clr=%.3f (d%+.3f)  ",e$denoiser[i],e$raw_auc[i],e$clr_auc[i],e$clr_auc[i]-e$raw_auc[i]))
    cat("\n") }
  rm(M); invisible(gc()) }
res<-do.call(rbind,rows); res$delta<-res$clr_auc-res$raw_auc
write.csv(res,file.path(out_dir,"generality_denoisers_results.csv"),row.names=FALSE)

#summary
agg<-res|>group_by(dataset,denoiser)|>summarise(raw=mean(raw_auc),clr=mean(clr_auc),
   delta=mean(delta),pos=sum(delta>0),n=n(),.groups="drop")|>
   mutate(dataset=factor(dataset,levels=c("Paul","PBMC","TM","Tusi")),
          denoiser=factor(denoiser,levels=c("raw","knn","svd","alra","magic")))|>arrange(dataset,denoiser)
cat("\n\n==== CLR vs raw readout, by denoiser (mean over",length(SEEDS),"seeds) ====\n")
print(as.data.frame(agg),row.names=FALSE,digits=3)
cat(sprintf("\nCLR > raw in %d / %d denoiser x dataset x seed cells\n",sum(res$delta>0),nrow(res)))

#figure
NAVY<-"#213C51";BLUE<-"#6594B1";LBLUE<-"#A9C5DA";ROSE<-"#C97FA8";PINK<-"#DDAED3"
theme_annice<-function(base_size=10){theme_minimal(base_size=base_size)+theme(
  panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),
  plot.background=element_blank(),axis.line=element_line(colour="black"),axis.ticks=element_line(colour="black"),
  axis.text=element_text(colour="black"),plot.title=element_text(hjust=0.5,face="bold"),
  strip.background=element_blank(),strip.text=element_text(face="bold"),legend.position="none")}
lab<-c(raw="no imp.",knn="kNN",svd="SVD",alra="ALRA",magic="MAGIC")
figd<-agg|>mutate(denoiser=factor(lab[as.character(denoiser)],levels=c("no imp.","kNN","SVD","ALRA","MAGIC")))
fig<-ggplot(figd,aes(denoiser,delta,fill=denoiser))+geom_col(width=.72)+
  geom_hline(yintercept=0,linewidth=.3)+facet_wrap(~dataset,nrow=1)+
  scale_fill_manual(values=c("no imp."=NAVY,"kNN"=BLUE,"SVD"=LBLUE,"ALRA"=ROSE,"MAGIC"=PINK))+
  labs(x="denoiser",y="CLR − raw  TF-target AUC gain",
       title="CLR null-calibration improves recovery under every denoiser (mean over 5 seeds)")+
  theme_annice(10)+theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(fig_dir,"fig_generality_denoisers.png"),fig,width=11,height=3.4,dpi=220,bg="white")
cat("\nwrote generality_denoisers_results.csv and fig_generality_denoisers.png\n")
