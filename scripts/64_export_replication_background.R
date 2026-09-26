# Export the sparse author ATAC matrix for genome-wide pseudobulk normalization.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
object <- readRDS(args[1])
counts <- slot(slot(object, "assays")[["ATAC"]], "counts")
rm(object)
gc()
writeLines(rownames(counts), file.path(out, "all_ATAC_features.txt"))
writeLines(colnames(counts), file.path(out, "all_ATAC_barcodes.txt"))
writeLines(as.character(dim(counts)), file.path(out, "all_ATAC_shape.txt"))
writeBin(counts@i, file.path(out, "all_ATAC_i.bin"), size = 4, endian = "little")
writeBin(counts@p, file.path(out, "all_ATAC_p.bin"), size = 4, endian = "little")
con <- file(file.path(out, "all_ATAC_x.bin"), "wb")
for (start in seq.int(1L, length(counts@x), by = 10000000L)) {
    idx <- start:min(start + 9999999L, length(counts@x))
    values <- counts@x[idx]
    stopifnot(all(is.finite(values)), all(values >= 0), all(values == floor(values)), max(values) < .Machine$integer.max)
    writeBin(as.integer(values), con, size = 4, endian = "little")
}
close(con)
cat("Exported", nrow(counts), "peaks;", ncol(counts), "nuclei;", length(counts@x), "nonzero entries\n")
