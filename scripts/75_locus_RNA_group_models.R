# Same matched nuclei and RNA library offsets used for the chr16 ATAC follow-up.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4)
root <- args[1]; work <- args[2]; locus <- args[3]; out <- args[4]
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(edgeR))
background <- as.matrix(readMM(file.path(work,"sensitivity_RNA_pseudobulk.mtx")))
rownames(background) <- readLines(file.path(work,"candidate_RNA_genes.txt"))
local <- as.matrix(readMM(file.path(locus,"locus_RNA_pseudobulk.mtx")))
rownames(local) <- readLines(file.path(locus,"locus_RNA_genes.txt"))
targets <- rownames(local)
counts <- rbind(background,local[!rownames(local) %in% rownames(background),,drop=FALSE])
cols <- read.delim(file.path(root,"results/replication_sensitivity/pseudobulk_columns.tsv"))
output <- list()
for(cfg in c("C169","C170","C171","C197","C198","C199")) {
 c <- cols[cols$config==cfg,];c$donor <- factor(c$donor);c$group <- factor(c$group,levels=c("low","high"))
 design <- model.matrix(~donor+group,c)
 y <- DGEList(counts[,c$column+1,drop=FALSE],lib.size=c$raw_RNA_size)
 bd <- sapply(levels(c$donor),function(d) rowSums(y$counts[,c$donor==d,drop=FALSE]))
 keep <- rowSums(y$counts)>=10 & rowSums(bd>0)>=2
 y <- y[keep,,keep.lib.sizes=TRUE]
 y <- estimateDisp(y,design,robust=TRUE)
 fit <- glmQLFit(y,design,robust=TRUE)
 t <- topTags(glmQLFTest(fit,coef=which(colnames(design)=="grouphigh")),n=Inf,sort.by="none")$table
 t$gene <- rownames(t); t <- t[t$gene %in% targets,];t$FC <- 2^t$logFC;t$config <- cfg;t$n_donors <- nlevels(c$donor)
 t$q_local_targets <- p.adjust(t$PValue,"BH");t$n_local_targets <- nrow(t)
 output[[cfg]] <- t
}
write.table(do.call(rbind,output),file.path(out,"RNA_high_low_models.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_model_session.txt"))
