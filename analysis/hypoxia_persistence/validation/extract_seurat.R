# Extract HIF-gene counts + metadata from a Seurat RDS without keeping the full object in Python.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: extract_seurat.R <input.rds> <out_prefix> [gene_file]")
}
infile <- args[[1]]
prefix <- args[[2]]
library(Seurat)
library(SeuratObject)
obj <- readRDS(infile)
DefaultAssay(obj) <- DefaultAssay(obj)
counts <- tryCatch(
  GetAssayData(obj, layer = "counts"),
  error = function(e) GetAssayData(obj, slot = "counts")
)
genes <- rownames(counts)
if (length(args) >= 3) {
  genes <- scan(args[[3]], what = character(), quiet = TRUE)
}
present <- intersect(genes, rownames(counts))
# also try TITLE case / aliases
if (length(present) < length(genes)) {
  rn <- rownames(counts)
  up <- toupper(rn)
  for (g in setdiff(genes, present)) {
    hit <- rn[up == toupper(g)]
    if (length(hit)) present <- unique(c(present, hit[[1]]))
  }
}
sub <- counts[present, , drop = FALSE]
meta <- slot(obj, "meta.data")
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
Matrix::writeMM(sub, paste0(prefix, "_counts.mtx"))
write.table(data.frame(gene = rownames(sub)), paste0(prefix, "_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(cell = colnames(sub)), paste0(prefix, "_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(meta, paste0(prefix, "_meta.csv"), row.names = TRUE)
cat("wrote", prefix, "genes", length(present), "cells", ncol(sub), "\n")
cat("meta columns:", paste(colnames(meta), collapse = ", "), "\n")
