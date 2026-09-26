# Donor-blocked QL models, with genome-wide ATAC normalization.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=3)
root <- args[1]; work <- args[2]; out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(Matrix))
configs <- read.delim(file.path(out,"configurations.tsv"),check.names=FALSE)
columns <- read.delim(file.path(out,"pseudobulk_columns.tsv"),check.names=FALSE)
features <- readLines(file.path(work,"all_ATAC_features.txt"))
nfeat <- length(features)
overlap <- read.delim(file.path(work,"candidate_peak_overlap.tsv"))
overlap <- subset(overlap,target_overlap_fraction>=.5 & author_overlap_fraction>=.5)
selected <- unique(overlap$author_peak)
idx <- match(selected,features)
stopifnot(!anyNA(idx))
annot <- read.delim(file.path(root,"results/enhancer_regions/peak_selectivity.tsv"))
fusion <- read.delim(file.path(root,"results/candidate_provenance/external_fusion_promoter_candidates.tsv"))
scope <- data.frame(peak=selected,all_candidates=TRUE,
    HSMM=selected %in% overlap$author_peak[overlap$target_peak %in% annot$peak],
    muscle_selective=selected %in% overlap$author_peak[overlap$target_peak %in% annot$peak[annot$n_other_strong==0]],
    fusion_promoter=selected %in% overlap$author_peak[overlap$target_peak %in% fusion$peak])
write.table(scope,file.path(out,"tested_peak_scopes.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
genes <- readLines(file.path(work,"candidate_RNA_genes.txt"))
rna <- as.matrix(readMM(file.path(work,"sensitivity_RNA_pseudobulk.mtx")))
rownames(rna) <- genes
stopifnot(ncol(rna)==nrow(columns))
bin <- file(file.path(work,"sensitivity_ATAC_pseudobulk.bin"),"rb")
atac_out <- gzfile(file.path(out,"ATAC_donor_models.tsv.gz"),"wt")
rna_out <- gzfile(file.path(out,"RNA_donor_models.tsv.gz"),"wt")
first_atac <- TRUE; first_rna <- TRUE; status <- list(); factors <- list()
requested <- if(length(args)>=4) strsplit(args[4],",")[[1]] else configs$config

fit_counts <- function(y,design,donors) {
    bydonor <- sapply(levels(donors),function(d) rowSums(y$counts[,donors==d,drop=FALSE]))
    keep <- rowSums(y$counts)>=10 & rowSums(bydonor>0)>=2
    y <- y[keep,,keep.lib.sizes=TRUE]
    if(nrow(y)<10) stop("Fewer than 10 count-supported features")
    y <- estimateDisp(y,design,robust=TRUE)
    fit <- glmQLFit(y,design,robust=TRUE)
    test <- glmQLFTest(fit,coef=which(colnames(design)=="grouphigh"))
    tab <- topTags(test,n=Inf,sort.by="none")$table
    list(table=tab,common_dispersion=y$common.dispersion)
}

for(cfg in requested) {
    c <- columns[columns$config==cfg,,drop=FALSE]
    nd <- length(unique(c$donor))
    if(nd<2) {
        status[[cfg]] <- data.frame(config=cfg,status="fewer_than_two_usable_donors",n_donors=nd,n_ATAC_tested=0)
        next
    }
    c$donor <- factor(c$donor); c$group <- factor(c$group,levels=c("low","high"))
    c$mean_log_RNA_centered <- c$mean_log_RNA-ave(c$mean_log_RNA,c$donor)
    design <- if(c$caliper[1]=="none") model.matrix(~donor+mean_log_RNA_centered+group,c) else model.matrix(~donor+group,c)
    if(qr(design)$rank<ncol(design) || nrow(design)<=ncol(design)) {
        status[[cfg]] <- data.frame(config=cfg,status="design_not_estimable",n_donors=nd,n_ATAC_tested=0)
        next
    }
    message(cfg," ",c$celltype[1]," ",c$time[1]," ",c$rule[1]," ",c$caliper[1]," N=",nd)
    stage <- "ATAC"; n_atac_written <- 0L
    tryCatch({
        counts <- vapply(c$column,function(j) {
            seek(bin,where=as.double(j)*nfeat*4,origin="start")
            readBin(bin,what="integer",n=nfeat,size=4,endian="little")
        },integer(nfeat))
        stopifnot(all(colSums(counts)==c$raw_ATAC_size))
        yfull <- DGEList(counts=counts)
        yfull <- calcNormFactors(yfull,method="TMM")
        fac <- c[,c("config","column","donor","group","n_nuclei")]
        fac$library_size <- yfull$samples$lib.size
        fac$TMM_factor <- yfull$samples$norm.factors
        factors[[cfg]] <- fac
        y <- yfull[idx,,keep.lib.sizes=TRUE]
        rownames(y$counts) <- selected
        result <- fit_counts(y,design,c$donor)
        tab <- result$table
        tab$peak <- rownames(tab);tab$config <- cfg;tab$n_donors <- nd
        tab$FC <- 2^tab$logFC
        sc <- scope[match(tab$peak,scope$peak),]
        for(name in names(scope)[-1]) {
            tab[[paste0("q_",name)]] <- NA_real_
            take <- sc[[name]]
            tab[[paste0("q_",name)]][take] <- p.adjust(tab$PValue[take],method="BH")
            tab[[paste0("n_",name)]] <- sum(take)
        }
        write.table(tab,atac_out,sep="\t",row.names=FALSE,quote=FALSE,col.names=first_atac,na="NA")
        first_atac <- FALSE
        n_atac_written <- nrow(tab); stage <- "RNA"
        yrna <- DGEList(counts=rna[,c$column+1,drop=FALSE],lib.size=c$raw_RNA_size)
        r <- fit_counts(yrna,design,c$donor)$table
        r$gene <- rownames(r);r$config <- cfg;r$n_donors <- nd;r$FC <- 2^r$logFC
        meaningful <- !r$gene %in% c("COQ8A","COQ9")
        r$q_myogenesis <- NA_real_
        r$q_myogenesis[meaningful] <- p.adjust(r$PValue[meaningful],method="BH")
        write.table(r,rna_out,sep="\t",row.names=FALSE,quote=FALSE,col.names=first_rna,na="NA")
        first_rna <- FALSE
        status[[cfg]] <- data.frame(config=cfg,status="ok",n_donors=nd,n_ATAC_tested=nrow(tab))
    },error=function(e) {
        status[[cfg]] <<- data.frame(config=cfg,status=paste0(stage,"_error: ",conditionMessage(e)),n_donors=nd,n_ATAC_tested=n_atac_written)
        message("ERROR ",cfg,": ",conditionMessage(e))
    })
    flush(atac_out);flush(rna_out)
    write.table(do.call(rbind,status),file.path(out,"model_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
}
close(bin);close(atac_out);close(rna_out)
write.table(do.call(rbind,status),file.path(out,"model_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,factors),file.path(out,"ATAC_normalization_factors.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_model_session.txt"))
cat("Completed",length(status),"configurations\n")
