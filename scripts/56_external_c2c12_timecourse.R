# Independent C2C12 differentiation reference: GSE224489.
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
input <- args[1]
out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
dat <- read.delim(gzfile(input), check.names = FALSE)
fields <- strsplit(dat$ID, "|", fixed = TRUE)
symbols <- vapply(fields, tail, character(1), n = 1)
mat <- as.matrix(dat[, -1])
rownames(mat) <- dat$ID
stopifnot(!anyNA(mat), all(mat >= 0), all(mat == round(mat)))
times <- vapply(colnames(mat), function(x) {
  if (grepl("_GM_GM_", x)) return("0")
  sub(".*_DM([0-9]+)h_[0-9]+$", "\\1", x)
}, character(1))
stopifnot(setequal(unique(times), c("0", "12", "24", "48", "60", "96")))
meta <- data.frame(sample = colnames(mat), time_h = factor(times, levels = c("0", "12", "24", "48", "60", "96")), row.names = colnames(mat))
stopifnot(all(table(meta$time_h) == 3))
write.table(meta, file.path(out, "sample_design.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
keep <- rowSums(mat) >= 10
dds <- DESeqDataSetFromMatrix(round(mat[keep, ]), meta, design = ~ time_h)
dds <- DESeq(dds, quiet = TRUE)
annotation <- data.frame(ID = dat$ID[keep], gene = symbols[keep])
normalized <- cbind(annotation, counts(dds, normalized = TRUE))
gz <- gzfile(file.path(out, "normalized_counts.tsv.gz"), "wt")
write.table(normalized, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)
results_all <- do.call(rbind, lapply(levels(meta$time_h)[-1], function(t) {
  r <- as.data.frame(results(dds, contrast = c("time_h", t, "0"), alpha = 0.05))
  r <- cbind(annotation, time_h = as.integer(t), r)
  r$FC <- 2^r$log2FoldChange
  r
}))
gz <- gzfile(file.path(out, "all_gene_time_contrasts.tsv.gz"), "wt")
write.table(results_all, gz, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
close(gz)
targets <- c("Coq8a", "Coq9", "Myod1", "Myog", "Myf5", "Cav3", "Csrp3", "Dmd", "Pfn2", "Rtn2", "Hspb8", "Crebrf", "Sgca", "Sgcd", "Lsp1", "Tnni2", "Tnnt3")
selected <- results_all[results_all$gene %in% targets, ]
write.table(selected, file.path(out, "candidate_time_contrasts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(data.frame(sample = names(sizeFactors(dds)), size_factor = sizeFactors(dds)), file.path(out, "size_factors.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), file.path(out, "R_session.txt"))
print(selected[selected$time_h %in% c(24, 48, 96), c("gene", "time_h", "FC", "pvalue", "padj")], row.names = FALSE)
