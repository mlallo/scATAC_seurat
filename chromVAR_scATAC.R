#!/usr/bin/env Rscript

# === Load arguments ===
args <- commandArgs(trailingOnly = TRUE)
atac_file <- args[1]
out_dir <- args[2]
file_prefix <- args[3]

# === Load libraries ===
library(Seurat)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(JASPAR2022)
library(TFBSTools)

register(MulticoreParam(7)) 
set.seed(2017)

# === Load Seurat object ===
cat("Reading Seurat object: ", atac_file, "\n")
atac_merged <- readRDS(atac_file)
DefaultAssay(atac_merged) <- "peaks"

# === Add patient.id to metadata ===
atac_merged$patient.id <- sub("^([^_]+)_.*", "\\1", rownames(atac_merged@meta.data))


# === List of unique patients ===
patients <- unique(atac_merged$patient.id)
cat("Unique patients: ", paste(patients, collapse = ", "), "\n")

# === Load motifs ===
opts <- list(species = 9606, all_versions = FALSE)
motifs <- getMatrixSet(JASPAR2022, opts)

# === Loop through each patient ===
for (pid in patients) {
  cat("Processing patient:", pid, "\n")
  
  # Subset Seurat object
  atac_patient <- subset(atac_merged, subset = patient.id == pid)
  
  # Extract count matrix and ranges
  peak_counts <- GetAssayData(atac_patient, assay = "peaks", layer = "counts")
  peak_ranges <- granges(atac_patient)
  
  # Create chromVAR input
  chromvar_input <- SummarizedExperiment(
    assays = list(counts = peak_counts),
    rowRanges = peak_ranges,
    colData = atac_patient@meta.data
  )
  
  # === Chromosome filtering ===
  standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
  keep_rows <- seqnames(rowRanges(chromvar_input)) %in% standard_chroms
  chromvar_input <- chromvar_input[keep_rows, ]
  rowRanges(chromvar_input) <- dropSeqlevels(
    rowRanges(chromvar_input),
    setdiff(seqlevels(rowRanges(chromvar_input)), standard_chroms),
    pruning.mode = "coarse"
  )
  
  # === Fix seqinfo and trim ===
  seqinfo(rowRanges(chromvar_input)) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  rowRanges(chromvar_input) <- trim(rowRanges(chromvar_input))
  
  # === GC Bias and depth ===
  chromvar_input <- addGCBias(chromvar_input, genome = BSgenome.Hsapiens.UCSC.hg19)
  colData(chromvar_input)$depth <- Matrix::colSums(assay(chromvar_input, "counts"))
  cat("Summary of depth for", pid, ":\n")
  print(summary(colData(chromvar_input)$depth))
  
  # === QC Plotting: Depth and In-Peak Fraction ===

depth_vec <- colData(chromvar_input)$depth
in_peaks_vec <- Matrix::colSums(assay(chromvar_input)) / depth_vec

cat("Cell count before filtering: ", length(depth_vec), "\n")
cat("In-peak fraction summary:\n")
print(summary(in_peaks_vec))

# Plot: depth distribution
depth_plot <- ggplot(data.frame(depth = depth_vec), aes(x = depth)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  scale_x_log10() +
  theme_minimal() +
  ggtitle(paste("Depth Distribution –", pid)) +
  xlab("Fragments per Cell (log10)") + ylab("Cell Count")

# Plot: in-peak fraction
peak_plot <- ggplot(data.frame(in_peaks = in_peaks_vec), aes(x = in_peaks)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black") +
  theme_minimal() +
  ggtitle(paste("In-Peak Fraction –", pid)) +
  xlab("Fraction of Reads in Peaks") + ylab("Cell Count")

# Save plots
qc_plot_dir <- file.path(out_dir, "qc_plots")
dir.create(qc_plot_dir, showWarnings = FALSE)

ggsave(file.path(qc_plot_dir, paste0(file_prefix, "_", pid, "_depth_hist.png")), plot = depth_plot, width = 6, height = 4)
ggsave(file.path(qc_plot_dir, paste0(file_prefix, "_", pid, "_inpeaks_hist.png")), plot = peak_plot, width = 6, height = 4)

  # === Filtering ===
  chromvar_filtered <- filterSamples(chromvar_input, min_depth = 1500, min_in_peaks = 0.15)
  chromvar_filtered <- filterPeaks(chromvar_filtered)
  chromvar_filtered <- chromvar_filtered[!is.na(rowData(chromvar_filtered)$bias), ]
  
  # === Match motifs and compute deviations ===
  motif_ix <- matchMotifs(motifs, chromvar_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)
  dev <- computeDeviations(object = chromvar_filtered, annotations = motif_ix)
  
  # === Save results per patient ===
  out_file <- file.path(out_dir, paste0(file_prefix, "_", pid, "_chromVAR.RData"))
  save(dev, motifs, motif_ix, chromvar_filtered, file = out_file)
  cat("Saved chromVAR for", pid, "→", out_file, "\n\n")
}
