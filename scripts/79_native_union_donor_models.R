# Donor-blocked models of every native 25-gene candidate interval.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=3)
root <- args[1]; work <- args[2]; out <- args[3]
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(Matrix))
base <- file.path(root,"results/replication_sensitivity")
configs <- read.delim(file.path(base,"configurations.tsv"),check.names=FALSE)
columns <- read.delim(file.path(base,"pseudobulk_columns.tsv"),check.names=FALSE)
features <- readLines(file.path(work,"all_ATAC_features.txt"))
nfeat <- length(features)
native <- read.delim(file.path(out,"native_peak_membership.tsv"))
background <- read.delim(file.path(base,"tested_peak_scopes.tsv"))$peak
selected <- unique(c(background,native$peak))
idx <- match(selected,features)
stopifnot(!anyNA(idx))
norm <- read.delim(file.path(base,"ATAC_normalization_factors.tsv"))
requested <- if(length(args)>=4) strsplit(args[4],",")[[1]] else configs$config
bin <- file(file.path(work,"sensitivity_ATAC_pseudobulk.bin"),"rb")
dest <- gzfile(file.path(out,"ATAC_native_models.tsv.gz"),"wt")
first <- TRUE; status <- list()
for(cfg in requested) {
 c <- columns[columns$config==cfg,,drop=FALSE]
 nd <- length(unique(c$donor))
 if(nd<2) {status[[cfg]] <- data.frame(config=cfg,status="fewer_than_two_usable_donors");next}
 c$donor <- factor(c$donor); c$group <- factor(c$group,levels=c("low","high"))
 c$mean_log_RNA_centered <- c$mean_log_RNA-ave(c$mean_log_RNA,c$donor)
 design <- if(c$caliper[1]=="none") model.matrix(~donor+mean_log_RNA_centered+group,c) else model.matrix(~donor+group,c)
 if(qr(design)$rank<ncol(design) || nrow(design)<=ncol(design)) {status[[cfg]] <- data.frame(config=cfg,status="design_not_estimable");next}
 tryCatch({
  counts <- vapply(c$column,function(j) {
   seek(bin,where=as.double(j)*nfeat*4,origin="start")
   readBin(bin,what="integer",n=nfeat,size=4,endian="little")[idx]
  },integer(length(idx)))
  fac <- norm$TMM_factor[match(c$column,norm$column)]
  if(anyNA(fac)) stop("Missing genome-wide TMM factor")
  y <- DGEList(counts=counts,lib.size=c$raw_ATAC_size,norm.factors=fac)
  rownames(y$counts) <- selected
  bydonor <- sapply(levels(c$donor),function(d) rowSums(y$counts[,c$donor==d,drop=FALSE]))
  keep <- rowSums(y$counts)>=10 & rowSums(bydonor>0)>=2
  y <- y[keep,,keep.lib.sizes=TRUE]
  if(nrow(y)<10) stop("Fewer than 10 count-supported background features")
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  test <- glmQLFTest(fit,coef=which(colnames(design)=="grouphigh"))
  tab <- topTags(test,n=Inf,sort.by="none")$table
  tab$peak <- rownames(tab)
  tab <- tab[tab$peak %in% native$peak,,drop=FALSE]
  tab$config <- cfg; tab$n_donors <- nd; tab$FC <- 2^tab$logFC
  tab$q_native_union <- p.adjust(tab$PValue,method="BH")
  tab$n_native_tested <- nrow(tab)
  write.table(tab,dest,sep="\t",row.names=FALSE,quote=FALSE,col.names=first,na="NA")
  first <- FALSE
  status[[cfg]] <- data.frame(config=cfg,status="ok")
  message(cfg," N=",nd," native peaks=",nrow(tab))
 },error=function(e) {status[[cfg]] <<- data.frame(config=cfg,status=paste0("error: ",conditionMessage(e)));message(cfg," ",conditionMessage(e))})
 flush(dest)
}
close(bin);close(dest)
write.table(do.call(rbind,status),file.path(out,"model_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,"R_model_session.txt"))
