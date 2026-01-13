#!/usr/bin/env Rscript

# Load required libraries
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(hdf5r)

# Define Base Directories for ATAC samples
base_dirs <- list(
  AML = "/atac/",
  CTRL = "/atac/"
)

out_base_dir <- "/atac_analysis"

# Create output directories if they don't exist
dir.create(out_base_dir, showWarnings = FALSE, recursive = TRUE)

# Output paths
output_dir <- file.path(out_base_dir, "objects")
figure_dir <- file.path(out_base_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


message("Starting scATAC-seq processing...")

# Get List of Sample Folders
samples <- list()
for (type in names(base_dirs)) {
  dirs <- list.dirs(base_dirs[[type]], recursive = FALSE, full.names = FALSE)
  if (type == "AML") {
    samples[[type]] <- dirs[grepl("^AML", dirs)]
  } else if (type == "CTRL") {
    samples[[type]] <- dirs[grepl("^CTRL", dirs)]
  }
}

# Initialize List to Store Seurat Objects
seurat_objects <- list()

# Process Each ATAC Sample
for (type in names(samples)) {
  for (sample in samples[[type]]) {
    message(paste("Processing sample:", sample))

    # Define paths to required files
    base_dir <- base_dirs[[type]]
    peak_matrix_file <- file.path(base_dir, sample, "outs", "filtered_peak_bc_matrix.h5")
    metadata_file <- file.path(base_dir, sample, "outs", "singlecell.csv")
    fragments_file <- file.path(base_dir, sample, "outs", "fragments.tsv.gz")

    # Ensure all required files exist
    if (!file.exists(peak_matrix_file) | !file.exists(metadata_file) | !file.exists(fragments_file)) {
      message(paste("Skipping", sample, "- Missing required files."))
      next
    }

    # Read ATAC Peak Matrix and Metadata
    message(paste("Reading peak matrix and metadata for:", sample))
    counts <- Read10X_h5(filename = peak_matrix_file)
    metadata <- read.csv(metadata_file, header = TRUE, row.names = 1)

    # Create Chromatin Assay
    message(paste("Creating Chromatin Assay for:", sample))
    chrom_assay <- CreateChromatinAssay(
      counts = counts, sep = c(":", "-"),
      fragments = fragments_file,
      min.cells = 10, min.features = 200
    )

    # Create Seurat Object
    seurat_objects[[sample]] <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)

    message(paste("Finished processing sample:", sample))
  }
}

# Extract Gene Annotations
message("Extracting gene annotations...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "grch38.86"

# Add Annotations & Compute QC Metrics
for (sample in names(seurat_objects)) {
  message(paste("Computing QC metrics for:", sample))
  Annotation(seurat_objects[[sample]]) <- annotations
  seurat_objects[[sample]] <- NucleosomeSignal(seurat_objects[[sample]])
  seurat_objects[[sample]] <- TSSEnrichment(seurat_objects[[sample]], fast = FALSE)
  seurat_objects[[sample]]$pct_reads_in_peaks <- seurat_objects[[sample]]$peak_region_fragments / seurat_objects[[sample]]$passed_filters * 100
  seurat_objects[[sample]]$blacklist_ratio <- seurat_objects[[sample]]$blacklist_region_fragments / seurat_objects[[sample]]$peak_region_fragments
}

# Save Seurat Objects after QC
saveRDS(seurat_objects, file.path(output_dir, "seurat_objects_qc.rds"))
message("Saved Seurat objects after QC.")

# Load scRNA-seq Reference Data
message("Loading scRNA-seq reference data...")
ubtf_rna <- readRDS("/seurat_obj_aml.rds")

# Quality Control Filtering
filtered_objects <- list()
for (sample in names(seurat_objects)) {
  message(paste("Filtering low-quality cells for:", sample))
  filtered_objects[[sample]] <- subset(seurat_objects[[sample]],
    subset = pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
             nucleosome_signal < 4 & TSS.enrichment > 3
  )
}

# Save QC-filtered objects
saveRDS(filtered_objects, file.path(output_dir, "seurat_objects_filtered_qc.rds"))
message("Saved QC-filtered Seurat objects.")

# Normalization & Dimensional Reduction
for (sample in names(filtered_objects)) {
  message(paste("Running TF-IDF normalization for:", sample))
  filtered_objects[[sample]] <- RunTFIDF(filtered_objects[[sample]])
  filtered_objects[[sample]] <- FindTopFeatures(filtered_objects[[sample]], min.cutoff = 'q0')
  filtered_objects[[sample]] <- RunSVD(filtered_objects[[sample]])
}

# Clustering & UMAP
for (sample in names(filtered_objects)) {
  message(paste("Running UMAP and clustering for:", sample))
  filtered_objects[[sample]] <- RunUMAP(filtered_objects[[sample]], reduction = 'lsi', dims = 2:30)
  filtered_objects[[sample]] <- FindNeighbors(filtered_objects[[sample]], reduction = 'lsi', dims = 2:30)
  filtered_objects[[sample]] <- FindClusters(filtered_objects[[sample]], verbose = FALSE, algorithm = 3)
}

# Save after clustering
saveRDS(filtered_objects, file.path(output_dir, "seurat_objects_after_clustering.rds"))
message("Saved Seurat objects after clustering.")

# Function to save clustering UMAP plots
save_clustering_plot <- function(seurat_obj, sample_name, output_folder) {
  plot_file <- file.path(output_folder, paste0(sample_name, "_clustering.png"))
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
    ggtitle(paste("UMAP Clustering:", sample_name)) +
    theme_minimal()
  
  ggsave(plot_file, plot = p, width = 6, height = 5, dpi = 300)
}

# Save clustering figures for each sample
for (sample in names(filtered_objects)) {
  message(paste("Saving UMAP clustering plot for:", sample))
  save_clustering_plot(filtered_objects[[sample]], sample, figure_dir)
}

message("All clustering figures saved.")

## has an error for the below

# Ensure scRNA-seq uses the RNA assay
DefaultAssay(ubtf_rna) <- "RNA"

# Initialize list to store labeled scATAC objects
filtered_labeled_objects <- list()

# Cell Type Labeling using scRNA Reference
for (sample in names(filtered_objects)) {
  message(paste("Labeling cell types for:", sample))

  # Extract scATAC-seq sample
  atac_sample <- filtered_objects[[sample]]

  # Check if gene activity exists, if not, create it
  if (!"RNA" %in% Assays(atac_sample)) {
    message(paste("Creating GeneActivity matrix for:", sample))
    atac_sample[["RNA"]] <- CreateAssayObject(counts = GeneActivity(atac_sample))
  }

  # Set RNA as the active assay for scATAC-seq
  DefaultAssay(atac_sample) <- "RNA"

  # Find shared features (common genes between RNA & ATAC)
  common_features <- intersect(rownames(ubtf_rna), rownames(atac_sample[["RNA"]]))
  if (length(common_features) == 0) {
    message(paste("Skipping", sample, "- No shared features found between RNA and ATAC datasets!"))
    next
  }

  message(paste("Shared features found for", sample, ":", length(common_features)))

  # Find transfer anchors using only shared features
  transfer.anchors <- FindTransferAnchors(
    reference = ubtf_rna,
    query = atac_sample,
    features = common_features,  # Use only shared features
    reduction = "cca",
    dims = 1:30
  )

  # Transfer RNA-seq labels to scATAC-seq cells
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = ubtf_rna$broad_cell_type,
    weight.reduction = atac_sample[['lsi']],  # Use LSI reduction from scATAC-seq
    dims = 2:30
  )

  # Add predicted labels as metadata
  atac_sample <- AddMetaData(object = atac_sample, metadata = predicted.labels)

  # Store labeled object
  filtered_labeled_objects[[sample]] <- atac_sample
} 

# Save Final Processed Seurat Objects
saveRDS(filtered_labeled_objects, file.path(output_dir, "scATAC_processed.rds"))
message("Final processed scATAC-seq objects saved.")

message("Processing complete!")
