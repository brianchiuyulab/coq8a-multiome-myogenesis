# Inspect the author Seurat object without replacing cell-type annotations.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
input <- args[1]
root <- args[2]
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(GenomicRanges))
cat("Reading author object:", input, "\n")
flush.console()
object <- readRDS(input)
cat("Loaded class:", class(object), "\n")
flush.console()
metadata <- slot(object, "meta.data")
metadata$barcode <- rownames(metadata)
write.table(metadata, file.path(out, "author_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
assays <- slot(object, "assays")
inventory <- lapply(names(assays), function(name) {
  assay <- assays[[name]]
  counts <- if ("counts" %in% names(attributes(assay))) slot(assay, "counts") else NULL
  data.frame(assay = name, class = class(assay)[1], n_features = if (is.null(counts)) NA else nrow(counts),
             n_nuclei = if (is.null(counts)) NA else ncol(counts),
             counts_class = if (is.null(counts)) NA else class(counts)[1])
})
write.table(do.call(rbind, inventory), file.path(out, "assay_inventory.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
capture.output(lapply(assays, function(x) names(attributes(x))), file = file.path(out, "assay_slots.txt"))
print(do.call(rbind, inventory))
flush.console()
stopifnot("RNA" %in% names(assays))
rna <- slot(assays[["RNA"]], "counts")
writeLines(rownames(rna), file.path(out, "author_RNA_features.txt"))
print(head(rownames(rna)))
print(rownames(rna)[grepl("COQ8|ADCK3", rownames(rna))])
stopifnot("COQ8A" %in% rownames(rna))
atac_names <- names(assays)[vapply(assays, function(x) identical(class(x)[1], "ChromatinAssay"), logical(1))]
stopifnot(length(atac_names) == 1)
atac <- slot(assays[[atac_names]], "counts")
stopifnot(!anyDuplicated(colnames(rna)), !anyDuplicated(colnames(atac)),
          setequal(colnames(rna), colnames(atac)), setequal(colnames(rna), rownames(metadata)))
atac <- atac[, colnames(rna), drop = FALSE]
metadata <- metadata[colnames(rna), , drop = FALSE]
position <- slot(assays[[atac_names]], "positionEnrichment")
capture.output(lapply(position, function(x) list(class = class(x), dim = dim(x))),
               file = file.path(out, "position_enrichment_inventory.txt"))
metadata$COQ8A_raw_UMI <- as.numeric(rna["COQ8A", ])
metadata$CAV3_raw_UMI <- if ("CAV3" %in% rownames(rna)) as.numeric(rna["CAV3", ]) else NA
metadata$raw_RNA_library_size <- Matrix::colSums(rna)
metadata$raw_ATAC_library_size <- Matrix::colSums(atac)
write.table(metadata, file.path(out, "nucleus_inventory.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(rownames(atac), file.path(out, "author_ATAC_features.txt"))
writeLines(rownames(rna), file.path(out, "author_RNA_features.txt"))
gene_space <- read.delim(file.path(root, "results/tables/gene_search_space.tsv"))$gene
genes <- intersect(unique(c("COQ8A", "COQ9", gene_space)), rownames(rna))
saveRDS(rna[genes, , drop = FALSE], file.path(out, "candidate_RNA_counts.rds"))
Matrix::writeMM(rna[genes, , drop = FALSE], file.path(out, "candidate_RNA_counts.mtx"))
writeLines(genes, file.path(out, "candidate_RNA_genes.txt"))
writeLines(colnames(rna), file.path(out, "barcodes.txt"))
old <- read.delim(file.path(root, "results/unstratified/ATAC_all_5097.tsv"))
parse_peaks <- function(ids) {
  pieces <- strsplit(ids, "[:-]")
  stopifnot(all(lengths(pieces) == 3))
  GRanges(vapply(pieces, `[`, character(1), 1),
          IRanges(as.integer(vapply(pieces, `[`, character(1), 2)) + 1L,
                  as.integer(vapply(pieces, `[`, character(1), 3))))
}
target <- parse_peaks(old$peak)
author <- parse_peaks(rownames(atac))
hits <- findOverlaps(target, author, ignore.strand = TRUE)
overlap_bp <- width(pintersect(target[queryHits(hits)], author[subjectHits(hits)], ignore.strand = TRUE))
mapping <- data.frame(target_peak = old$peak[queryHits(hits)], author_peak = rownames(atac)[subjectHits(hits)],
                      overlap_bp = overlap_bp, target_width = width(target)[queryHits(hits)],
                      author_width = width(author)[subjectHits(hits)])
mapping$target_overlap_fraction <- mapping$overlap_bp / mapping$target_width
mapping$author_overlap_fraction <- mapping$overlap_bp / mapping$author_width
write.table(mapping, file.path(out, "candidate_peak_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
idx <- sort(unique(subjectHits(hits)))
saveRDS(atac[idx, , drop = FALSE], file.path(out, "candidate_ATAC_counts.rds"))
Matrix::writeMM(atac[idx, , drop = FALSE], file.path(out, "candidate_ATAC_counts.mtx"))
writeLines(rownames(atac)[idx], file.path(out, "candidate_ATAC_peaks.txt"))
capture.output(sessionInfo(), file = file.path(out, "R_session.txt"))
cat("Export complete. Nuclei:", nrow(metadata), "COQ8A >=3:", sum(metadata$COQ8A_raw_UMI >= 3), "\n")
