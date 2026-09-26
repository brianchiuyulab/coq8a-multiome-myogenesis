# Sensitivity of selected candidate effects to covariates and donor composition.
args <- commandArgs(trailingOnly=TRUE)
options(warn=1)
stopifnot(length(args)==3)
root <- args[1]; work <- args[2]; out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(edgeR))
base <- file.path(root,"results/replication_sensitivity")
columns <- read.delim(file.path(base,"pseudobulk_columns.tsv"))
nf <- read.delim(file.path(base,"ATAC_normalization_factors.tsv"))
features <- readLines(file.path(work,"all_ATAC_features.txt"))
scope <- read.delim(file.path(base,"tested_peak_scopes.tsv"))
idx <- match(scope$peak,features)
focal <- c("chr3-8733222-8734203","chr1-211133720-211134678","chr16-1311478-1312392")
configs <- c(sprintf("C%03d",296:303),"C324","C328","C197","C198","C199","C169")
bin <- file(file.path(work,"sensitivity_ATAC_pseudobulk.bin"),"rb")
rows <- list(); donors_out <- list(); k <- 0L
for(cfg in configs) {
 c <- columns[columns$config==cfg,]; c <- merge(c,nf[,c("column","library_size","TMM_factor")],by="column",sort=FALSE)
 counts <- vapply(c$column,function(j){seek(bin,as.double(j)*length(features)*4,"start");readBin(bin,"integer",n=length(features),size=4,endian="little")[idx]},integer(length(idx)))
 rownames(counts) <- scope$peak
 for(peak in focal) {
  d <- c;d$peak <- peak;d$count <- counts[peak,];d$CPM <- d$count/(d$library_size*d$TMM_factor)*1e6
  donors_out[[length(donors_out)+1]] <- d
 }
 subsets <- list(all=seq_len(nrow(c)),exercise=which(c$donor %in% c("E","G","I","J")),rest=which(c$donor %in% c("L","N")))
 if(cfg %in% c("C296","C297","C198")) for(d in unique(c$donor)) subsets[[paste0("without_",d)]] <- which(c$donor!=d)
 for(sub in names(subsets)) for(covar in c("donor_only","donor_plus_RNA")) {
  ii <- subsets[[sub]];cc <- c[ii,];cc$donor <- factor(cc$donor);cc$group <- factor(cc$group,levels=c("low","high"))
  cc$RNA <- cc$mean_log_RNA-ave(cc$mean_log_RNA,cc$donor)
  nd <- nlevels(cc$donor)
  if(nd<2) next
  design <- if(covar=="donor_only") model.matrix(~donor+group,cc) else model.matrix(~donor+RNA+group,cc)
  if(qr(design)$rank<ncol(design) || nrow(design)<=ncol(design)) next
  tryCatch({
   factor_scale <- exp(mean(log(cc$TMM_factor)))
   y <- DGEList(counts[,ii,drop=FALSE],lib.size=cc$library_size*factor_scale,norm.factors=cc$TMM_factor/factor_scale)
   bd <- sapply(levels(cc$donor),function(d) rowSums(y$counts[,cc$donor==d,drop=FALSE]))
   y <- y[rowSums(y$counts)>=10 & rowSums(bd>0)>=2,,keep.lib.sizes=TRUE]
   y <- estimateDisp(y,design,robust=TRUE);fit <- glmQLFit(y,design,robust=TRUE)
   t <- topTags(glmQLFTest(fit,coef=which(colnames(design)=="grouphigh")),n=Inf,sort.by="none")$table
   sc <- scope[match(rownames(t),scope$peak),]
   t$q_fusion_promoter <- NA_real_;t$q_fusion_promoter[sc$fusion_promoter] <- p.adjust(t$PValue[sc$fusion_promoter],"BH")
   t$n_fusion_promoter <- sum(sc$fusion_promoter)
   t$peak <- rownames(t);t <- t[t$peak %in% focal,];t$config <- cfg;t$subset <- sub;t$model <- covar;t$n_donors <- nd;t$FC <- 2^t$logFC;t$design_condition_number <- kappa(design)
   k <- k+1;rows[[k]] <- t
  },error=function(e) message(cfg," ",sub," ",covar,": ",conditionMessage(e)))
 }
 message(cfg," complete")
}
close(bin)
write.table(do.call(rbind,rows),file.path(out,"candidate_model_stress_tests.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,donors_out),file.path(out,"candidate_counts_by_donor.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_stress_session.txt"))
