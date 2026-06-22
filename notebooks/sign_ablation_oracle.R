#sign-control for Table 5: |rho|, max(rho,0), CLR scoring, oracle over t=0:10, 5 seeds, random negatives
suppressPackageStartupMessages({library(Matrix); library(cvMAGIC); library(data.table); library(dplyr)})
data_dir <- normalizePath(file.path("..","data"), mustWork=TRUE)
out_dir  <- normalizePath(file.path("..","notebooks"), mustWork=FALSE)
NPCA<-20L; K<-30L; SEEDS<-101:105; T_GRID<-0:10

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

load_paul<-function(){cf<-fread(file.path(data_dir,"GSE72857_umitab.txt"),sep="\t",header=TRUE,data.table=FALSE);gn<-cf[[1]];cf<-as.matrix(cf[,-1]);rownames(cf)<-gn;M<-t(cf);storage.mode(M)<-"integer";M}
load_pbmc<-function(){d<-file.path(data_dir,"filtered_gene_bc_matrices","hg19");m<-readMM(file.path(d,"matrix.mtx"));gt<-read.table(file.path(d,"genes.tsv"),sep="\t");rownames(m)<-make.unique(gt$V2);m<-as.matrix(t(m));storage.mode(m)<-"integer";m}
load_tm<-function(){d<-file.path(data_dir,"droplet","Marrow-10X_P7_2");m<-readMM(file.path(d,"matrix.mtx"));gt<-read.table(file.path(d,"genes.tsv"),sep="\t");rownames(m)<-make.unique(gt$V1);m<-as.matrix(t(m));storage.mode(m)<-"integer";m}
load_tusi<-function(){csv<-file.path(data_dir,"GSM2388072_basal_bone_marrow.raw_umifm_counts.csv");df<-fread(csv,sep=",",header=TRUE,data.table=FALSE);meta<-c("cell_id","barcode","library_id","seq_run_id","pass_filter");M<-as.matrix(df[,setdiff(colnames(df),meta)]);rownames(M)<-df$cell_id;storage.mode(M)<-"integer";M}
specs<-list(
  Paul=list(load=load_paul,tf=tf_mouse,min_lib=500L,gene_pos=0.05),
  PBMC=list(load=load_pbmc,tf=tf_human,min_lib=500L,gene_pos=0.03),
  TM  =list(load=load_tm,  tf=tf_mouse,min_lib=500L,gene_pos=0.03),
  Tusi=list(load=load_tusi,tf=tf_tusi, min_lib=1000L,gene_pos=0.05))

build_panel<-function(M,tf,min_lib,gene_pos,seed,n_cells=1500L,n_genes=2000L){
  ls<-rowSums(M);keep<-which(ls>=min_lib&ls<=quantile(ls,0.99));set.seed(seed)
  ci<-sample(keep,min(n_cells,length(keep)));cn<-M[ci,,drop=FALSE]
  mk<-intersect(unique(c(names(tf),unlist(tf))),colnames(cn))
  gp<-colMeans(cn>0);cn<-cn[,gp>=gene_pos|colnames(cn)%in%mk,drop=FALSE]
  gv<-apply(cn,2,var);hv<-order(gv,decreasing=TRUE)[seq_len(min(n_genes,length(gv)))]
  X<-cn[,union(colnames(cn)[hv],mk),drop=FALSE];storage.mode(X)<-"integer";X}
norm_sqrt<-function(X){s<-rowSums(X);s[s==0]<-1;Xn<-X/s*mean(s);sqrt(Xn)}  #match residual_magic
make_pairs<-function(tf,gs,seed=7L){set.seed(seed);rows<-list()
  for(f in names(tf)){if(!(f%in%gs))next;tg<-intersect(tf[[f]],gs);if(!length(tg))next
    pool<-setdiff(gs,c(f,tg));ng<-sample(pool,min(length(pool),5L*length(tg)))
    rows[[f]]<-rbind(data.frame(tf=f,gene=tg,label=1L),data.frame(tf=f,gene=ng,label=0L))}
  do.call(rbind,rows)}
auc1<-function(scores,labels){r<-rank(scores);np<-sum(labels==1);nn<-sum(labels==0)
  if(np==0||nn==0)return(NA_real_);(sum(r[labels==1])-np*(np+1)/2)/(np*nn)}
val<-function(C,pr) mapply(function(a,b)C[a,b],pr$tf,pr$gene)
sc_absraw<-function(C,pr) abs(val(C,pr))                       #|rho|
sc_pos   <-function(C,pr) pmax(val(C,pr),0)                    #max(rho,0)
sc_clr1  <-function(C,pr){ mu<-rowMeans(C); sg<-apply(C,1,sd); sg[sg==0]<-1
  zi<-pmax(0,(val(C,pr)-mu[pr$tf])/sg[pr$tf]); zj<-pmax(0,(val(C,pr)-mu[pr$gene])/sg[pr$gene]); sqrt(zi^2+zj^2) }
variants<-list(`|rho|`=sc_absraw, `max(rho,0)`=sc_pos, CLR=sc_clr1)

rows<-list()
for(nm in names(specs)){ sp<-specs[[nm]]; M<-sp$load(); cat("===",nm,"===\n")
  for(s in SEEDS){ X<-build_panel(M,sp$tf,sp$min_lib,sp$gene_pos,seed=s); gs<-colnames(X); pr<-make_pairs(sp$tf,gs)
    Xn<-norm_sqrt(X); g<-magic_graph(Xn,npca=NPCA,k=K,ka=K); Y<-Xn; prev<-0L
    best<-setNames(rep(-Inf,length(variants)),names(variants))
    for(t in T_GRID){ while(prev<t){Y<-as.matrix(g$M%*%Y);prev<-prev+1L}
      imp<-if(t==0L) as.matrix(Xn) else Y; colnames(imp)<-gs
      C<-suppressWarnings(cor(imp,method="spearman")); C[is.na(C)]<-0
      for(v in names(variants)){ a<-auc1(variants[[v]](C,pr),pr$label); if(!is.na(a)&&a>best[v]) best[v]<-a } }
    rows[[paste(nm,s)]]<-data.frame(dataset=nm,seed=s,abs_rho=best[["|rho|"]],pos_rho=best[["max(rho,0)"]],clr=best[["CLR"]])
    cat(sprintf("  seed %d | |rho|=%.3f max(rho,0)=%.3f CLR=%.3f\n",s,best[["|rho|"]],best[["max(rho,0)"]],best[["CLR"]])) }
  rm(M); invisible(gc()) }
res<-do.call(rbind,rows); write.csv(res,file.path(out_dir,"sign_ablation_oracle_results.csv"),row.names=FALSE)

agg<-res|>group_by(dataset)|>summarise(abs_rho=mean(abs_rho),pos_rho=mean(pos_rho),clr=mean(clr),.groups="drop")|>
  mutate(d_sign=pos_rho-abs_rho,d_calib=clr-pos_rho,dataset=factor(dataset,levels=c("Paul","PBMC","TM","Tusi")))|>arrange(dataset)
cat("\n==== Table 5 (oracle over t=0:10, mean over 5 seeds, random negatives) ====\n")
print(as.data.frame(agg),row.names=FALSE,digits=3)
cat("\ncalibration adds positive gain in all seeds:",
    sum(res$clr>res$pos_rho),"/",nrow(res),"\n")
