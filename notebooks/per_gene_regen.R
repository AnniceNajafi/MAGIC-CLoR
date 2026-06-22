#per-gene Pearson adaptive t vs global MCV on Paul and PBMC: TF-target AUC (raw |Spearman|) and saturation fraction, 5 seeds
suppressPackageStartupMessages({library(Matrix); library(cvMAGIC); library(data.table); library(dplyr)})
data_dir <- normalizePath(file.path("..","data"), mustWork=TRUE)
out_dir  <- normalizePath(file.path("..","notebooks"), mustWork=FALSE)
NPCA<-20L; K<-30L; SEEDS<-101:105; TV<-1:10

tf_mouse<-list(
  Gata1=c("Klf1","Hba-a1","Hba-a2","Hbb-b1","Hbb-b2","Aqp1","Slc4a1","Gypa","Epor"),
  Gata2=c("Gata1","Klf1","Itga2b"),Klf1=c("Hba-a1","Hba-a2","Hbb-b1","Hbb-b2","Aqp1","Slc4a1"),
  Sfpi1=c("Cebpa","Mpo","Csf1r","Csf2ra","Csf2rb"),Cebpa=c("Mpo","Elane","Prtn3","Ctsg","Csf3r"),
  Cebpe=c("Mpo","Elane","Prtn3","Ctsg"),Gfi1=c("Csf1r","Elane"),Irf8=c("Cxcr4"),Tal1=c("Gata1","Klf1","Gata2"))
tf_human<-list(
  PAX5=c("CD19","MS4A1","CD79A","CD79B","BLNK","BANK1","BCL11A"),TBX21=c("IFNG","GZMB","CXCR3","GNLY","PRF1","NKG7"),
  EOMES=c("IFNG","GZMA","PRF1","GNLY","NKG7","KLRD1"),FOXP3=c("IL2RA","CTLA4","TNFRSF18"),
  SPI1=c("CD14","LYZ","CSF1R","FCGR3A","CST3","S100A8","S100A9"),IRF8=c("CCR7","FLT3","CST3","HLA-DRA"),
  CEBPB=c("CD14","LYZ","S100A8","S100A9"),LEF1=c("SELL","CCR7","TCF7"),KLF2=c("SELL","S1PR1","CCR7"),ID3=c("SELL","CCR7"))

load_paul<-function(){cf<-fread(file.path(data_dir,"GSE72857_umitab.txt"),sep="\t",header=TRUE,data.table=FALSE);gn<-cf[[1]];cf<-as.matrix(cf[,-1]);rownames(cf)<-gn;M<-t(cf);storage.mode(M)<-"integer";M}
load_pbmc<-function(){d<-file.path(data_dir,"filtered_gene_bc_matrices","hg19");m<-readMM(file.path(d,"matrix.mtx"));gt<-read.table(file.path(d,"genes.tsv"),sep="\t");rownames(m)<-make.unique(gt$V2);m<-as.matrix(t(m));storage.mode(m)<-"integer";m}
specs<-list(
  Paul=list(load=load_paul,tf=tf_mouse,min_lib=500L,gene_pos=0.05),
  PBMC=list(load=load_pbmc,tf=tf_human,min_lib=500L,gene_pos=0.03))

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
score_raw_auc<-function(imp,pr){ colnames(imp)<-colnames(imp); C<-suppressWarnings(cor(imp,method="spearman")); C[is.na(C)]<-0
  auc1(abs(mapply(function(a,b)C[a,b],pr$tf,pr$gene)),pr$label) }

rows<-list()
for(nm in names(specs)){ sp<-specs[[nm]]; M<-sp$load(); cat("===",nm,"===\n")
  for(s in SEEDS){ X<-build_panel(M,sp$tf,sp$min_lib,sp$gene_pos,seed=s); gs<-colnames(X); pr<-make_pairs(sp$tf,gs)
    g<-magic_graph(X,npca=NPCA,k=K)
    sel<-mcv_select_t(X,t_values=TV,npca=NPCA,k=K,seed=1L)
    imp_mcv<-magic_impute(X,t=sel$t,graph=g); colnames(imp_mcv)<-gs
    pg<-magic_per_gene(X,t_values=TV,loss="pearson",npca=NPCA,k=K,seed=1L)
    imp_pg<-pg$imputed; colnames(imp_pg)<-gs
    sat<-mean(pg$t_per_gene==max(TV))                 #fraction of genes saturating at t=max
    auc_mcv<-score_raw_auc(imp_mcv,pr); auc_pg<-score_raw_auc(imp_pg,pr)
    rows[[paste(nm,s)]]<-data.frame(dataset=nm,seed=s,t_global=sel$t,
      auc_global_mcv=auc_mcv,auc_per_gene_pearson=auc_pg,frac_saturated=sat)
    cat(sprintf("  seed %d | global MCV t=%d auc=%.3f | per-gene pearson auc=%.3f | sat=%.1f%%\n",
      s,sel$t,auc_mcv,auc_pg,100*sat)) }
  rm(M); invisible(gc()) }
res<-do.call(rbind,rows); write.csv(res,file.path(out_dir,"per_gene_regen_results.csv"),row.names=FALSE)

agg<-res|>group_by(dataset)|>summarise(global_mcv=mean(auc_global_mcv),per_gene_pearson=mean(auc_per_gene_pearson),
  frac_saturated=mean(frac_saturated),.groups="drop")
cat("\n==== Per-gene Pearson vs global MCV (mean over 5 seeds) ====\n")
print(as.data.frame(agg),row.names=FALSE,digits=3)
