# Export all measured genes in the unrestricted chr16 transcript-TSS neighborhood.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
root <- args[2]; out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
near <- read.delim(file.path(root,"results/candidate_identity/nearest_TSS_per_gene.tsv"))
near <- near[near$peak=="chr16:1311478-1312392" & near$midpoint_distance_bp<=100000,]
object <- readRDS(args[1])
rna <- slot(slot(object,"assays")[["RNA"]],"counts")
near$RNA_feature_available <- near$gene %in% rownames(rna)
genes <- unique(near$gene[near$RNA_feature_available])
writeMM(rna[genes,,drop=FALSE],file.path(out,"locus_RNA.mtx"))
writeLines(genes,file.path(out,"locus_RNA_genes.txt"))
writeLines(colnames(rna),file.path(out,"locus_RNA_barcodes.txt"))
write.table(near,file.path(out,"locus_RNA_inventory.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_export_session.txt"))
cat("Exported",length(genes),"genes for",ncol(rna),"nuclei\n")
