#raw |Spearman| vs CLR scoring on a fixed imputed matrix, pooled and per-TF AUC
suppressPackageStartupMessages({library(Matrix); library(cvMAGIC); library(data.table)})
data_dir <- normalizePath(file.path("..","data"), mustWork=TRUE)
NPCA<-20L; K<-30L; SEEDS<-101:103; T_DIFF<-3L

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
make_pairs<-function(tf,gs,seed=7L){set.seed(seed);rows<-list()
  for(f in names(tf)){if(!(f%in%gs))next;tg<-intersect(tf[[f]],gs);if(!length(tg))next
    pool<-setdiff(gs,c(f,tg));ng<-sample(pool,min(length(pool),5L*length(tg)))
    rows[[f]]<-rbind(data.frame(tf=f,gene=tg,label=1L),data.frame(tf=f,gene=ng,label=0L))}
  do.call(rbind,rows)}
auc1<-function(scores,labels){r<-rank(scores);np<-sum(labels==1);nn<-sum(labels==0)
  if(np==0||nn==0)return(NA_real_);(sum(r[labels==1])-np*(np+1)/2)/(np*nn)}
macro<-function(scores,pr){ tapply(seq_len(nrow(pr)),pr$tf,function(ix) auc1(scores[ix],pr$label[ix])) |> mean(na.rm=TRUE) }

val<-function(C,pr) mapply(function(a,b)C[a,b],pr$tf,pr$gene)
sc_absraw<-function(C,pr) abs(val(C,pr))                                   #|rho|
sc_pos   <-function(C,pr) pmax(val(C,pr),0)                                #max(rho,0)
clr_score<-function(Cb,pr,onesided=TRUE){ mu<-rowMeans(Cb); sg<-apply(Cb,1,sd); sg[sg==0]<-1
  zi<-(Cb[cbind(pr$tf,pr$gene)]-mu[pr$tf])/sg[pr$tf]; zj<-(Cb[cbind(pr$tf,pr$gene)]-mu[pr$gene])/sg[pr$gene]
  if(onesided){zi<-pmax(zi,0);zj<-pmax(zj,0)}; sqrt(zi^2+zj^2) }
sc_clr1  <-function(C,pr) clr_score(C,pr,TRUE)            #one-sided CLR (deployed)
sc_clr2  <-function(C,pr) clr_score(C,pr,FALSE)           #two-sided CLR
sc_clrabs<-function(C,pr) clr_score(abs(C),pr,TRUE)       #CLR on |rho|

variants<-list(`|rho|`=sc_absraw, `max(rho,0)`=sc_pos, `CLR|rho|`=sc_clrabs,
               `CLR-2sided`=sc_clr2, `CLR-1sided(deployed)`=sc_clr1)

res<-list()
for(nm in names(specs)){ sp<-specs[[nm]]; M<-sp$load(); print(paste0("=== ", nm, " ==="))
  for(s in SEEDS){ X<-build_panel(M,sp$tf,sp$min_lib,sp$gene_pos,seed=s); gs<-colnames(X); pr<-make_pairs(sp$tf,gs)
    g<-magic_graph(X,npca=NPCA,k=K)
    for(tt in c(0L,T_DIFF)){
      imp<-if(tt==0L) as.matrix(X) else { Y<-as.matrix(X); for(i in seq_len(tt)) Y<-g$M%*%Y; as.matrix(Y) }
      colnames(imp)<-gs; C<-suppressWarnings(cor(imp,method="spearman")); C[is.na(C)]<-0
      for(v in names(variants)){ sco<-variants[[v]](C,pr)
        res[[length(res)+1]]<-data.frame(dataset=nm,seed=s,t=tt,scoring=v,
          pooled=auc1(sco,pr$label), macro=macro(sco,pr)) } } } }
R<-do.call(rbind,res)
agg<-aggregate(cbind(pooled,macro)~dataset+t+scoring,R,mean)
write.csv(agg,"sign_ablation_results.csv",row.names=FALSE)

print(paste0("mean over ", length(SEEDS), " seeds"))
for(tt in c(0L,T_DIFF)){ print(paste0("t = ", tt))
  for(nm in names(specs)){ sub<-agg[agg$dataset==nm&agg$t==tt,]; line<-nm
    for(v in names(variants)){ r<-sub[sub$scoring==v,]
      line<-paste0(line, "  ", v, " pool=", round(r$pooled,3), " mac=", round(r$macro,3)) }
    print(paste0(line)) } }
