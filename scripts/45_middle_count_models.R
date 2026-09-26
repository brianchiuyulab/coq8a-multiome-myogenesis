# Raw-count association screen and stratified SCENT-style bootstrap follow-up.
args <- commandArgs(trailingOnly=TRUE)
directory <- args[1]
source_file <- args[2]
B <- if(length(args)>=3) as.integer(args[3]) else 1000L
suppressPackageStartupMessages(library(Matrix))
source(source_file)
meta <- read.delim(file.path(directory,'metadata.tsv'),check.names=FALSE)
genes <- readLines(file.path(directory,'genes.txt'))
peaks <- readLines(file.path(directory,'peaks.txt'))
rna <- as(readMM(file.path(directory,'rna.mtx')),'CsparseMatrix')
atac <- as(readMM(file.path(directory,'atac.mtx')),'CsparseMatrix')
rownames(rna)<-genes; rownames(atac)<-peaks
candidates<-read.delim(file.path(directory,'candidate_pairs.tsv'))
stopifnot(ncol(rna)==nrow(meta),ncol(atac)==nrow(meta))
models<-c('technical','technical_coq','technical_state','technical_state_coq')
design<-function(m,model,library=FALSE){
 z<-cbind(log1p(m$total_rna_umi),log1p(m$total_open_peaks),m$percent_mito)
 if(grepl('coq',model))z<-cbind(z,log1p(m$COQ8A_umi))
 if(grepl('state',model))z<-cbind(z,m$myogenesis_score)
 z<-scale(z); z[!is.finite(z)]<-0
 if(library)z<-cbind(model.matrix(~factor(m$gsm)),z) else z<-cbind(1,z)
 z
}
fit<-function(y,x){
 g<-glm.fit(x,y,family=poisson(),control=glm.control(maxit=50))
 if(!g$converged || any(!is.finite(g$coefficients)))return(NULL)
 bread<-tryCatch(solve(crossprod(x,x*g$fitted.values)),error=function(e)NULL)
 if(is.null(bread))return(NULL)
 meat<-crossprod(x*(y-g$fitted.values))
 se<-sqrt(diag(bread%*%meat%*%bread))
 list(beta=g$coefficients[1],se=se[1],converged=g$converged)
}
rows<-list();k<-0L
for(gsm in unique(meta$gsm)){
 ids<-which(meta$gsm==gsm); m<-meta[ids,]
 designs<-lapply(models,function(z)design(m,z));names(designs)<-models
 for(j in seq_len(nrow(candidates))){
  y<-as.numeric(rna[candidates$gene[j],ids]);a<-as.numeric(atac[candidates$peak[j],ids]>0)
  if(mean(y>0)<=.05 || mean(a>0)<=.05 || sum(y>0)<20 || sum(a)>length(a)-20)next
  for(model in models){
   f<-fit(y,cbind(a,designs[[model]]));if(is.null(f))next
   k<-k+1L;rows[[k]]<-data.frame(gsm=gsm,peak=candidates$peak[j],gene=candidates$gene[j],model=model,n=length(ids),beta=f$beta,se_robust=f$se,p=2*pnorm(-abs(f$beta/f$se)))
  }
 }
 cat(gsm,'count models complete\n');flush.console()
}
tab<-do.call(rbind,rows);write.table(tab,file.path(directory,'count_models_by_library.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
combined<-list();k<-0L
for(g in split(tab,interaction(tab$peak,tab$gene,tab$model,drop=TRUE))){
 if(nrow(g)<3)next
 w<-1/g$se_robust^2;b<-weighted.mean(g$beta,w);se<-sqrt(1/sum(w))
 k<-k+1L;combined[[k]]<-data.frame(peak=g$peak[1],gene=g$gene[1],model=g$model[1],n_libraries=nrow(g),beta=b,RNA_ratio=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p=2*pnorm(-abs(b/se)),positive_libraries=sum(g$beta>0))
}
out<-do.call(rbind,combined);out$q<-ave(out$p,out$model,FUN=function(p)p.adjust(p,'BH'))
write.table(out,file.path(directory,'count_models_combined.tsv'),sep='\t',row.names=FALSE,quote=FALSE)

# Fixed bootstrap budget and library-stratified resampling differ from SCENT's
# adaptive default. basic_p is imported unchanged from the pinned source.
targets<-data.frame(gene=c('LDB3','TNNC1'),peak=c('chr10:86674717-86675570','chr3:52446715-52447619'))
bootrows<-list();k<-0L
set.seed(208248)
for(j in seq_len(nrow(targets))){
 gene<-targets$gene[j];peak<-targets$peak[j]
 good<-unique(tab$gsm[tab$gene==gene & tab$peak==peak])
 ids<-which(meta$gsm%in%good);m<-meta[ids,]; y<-as.numeric(rna[gene,ids]);a<-as.numeric(atac[peak,ids]>0)
 for(model in c('technical','technical_state_coq')){
  x<-cbind(a,design(m,model,library=TRUE));base<-glm.fit(x,y,family=poisson());betas<-rep(NA_real_,B)
  groups<-split(seq_along(ids),m$gsm)
  for(b in seq_len(B)){
   ix<-unlist(lapply(groups,function(v)sample(v,length(v),replace=TRUE)),use.names=FALSE)
   z<-glm.fit(x[ix,,drop=FALSE],y[ix],family=poisson(),start=base$coefficients)
   if(z$converged)betas[b]<-z$coefficients[1]
   if(b%%250==0){cat(gene,model,b,'/',B,'bootstraps\n');flush.console()}
  }
  valid<-betas[is.finite(betas)];ci<-quantile(valid,c(.025,.975));k<-k+1L
  bootrows[[k]]<-data.frame(gene=gene,peak=peak,model=model,n_nuclei=length(ids),n_libraries=length(good),beta=base$coefficients[1],RNA_ratio=exp(base$coefficients[1]),CI_low=exp(ci[1]),CI_high=exp(ci[2]),p_boot=basic_p(base$coefficients[1],valid),B=length(valid))
  write.table(do.call(rbind,bootrows),file.path(directory,'scent_style_bootstrap_followup.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 }
}
