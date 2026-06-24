#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

source("github_utils.R")
cfg <- load_config()
set.seed(cfg$random_seed)

output_dir <- make_dir(file.path(cfg$result_dir, "01_clean_major_object"))
figure_dir <- make_dir(file.path(output_dir, "figures"))

message("Loading major-cell object: ", cfg$input_major_rds)
obj <- read_seurat(cfg$input_major_rds)
obj <- repair_legacy_spatial_slots(obj)

required_meta <- c(
  cfg$major_cluster_column, cfg$sample_column, cfg$age_column,
  cfg$position_column, cfg$x_column, cfg$y_column
)
require_columns(obj@meta.data, required_meta, "Seurat metadata")

cluster_values <- as.character(obj@meta.data[[cfg$major_cluster_column]])
keep_cells <- cluster_values != cfg$bad_cluster_label & !is.na(cluster_values)
message("Keeping ", sum(keep_cells), " cells; removing ", sum(!keep_cells), " cells.")

obj <- subset(obj, cells = colnames(obj)[keep_cells])
Idents(obj) <- cfg$major_cluster_column

if ("pca" %in% names(obj@reductions)) {
  dims <- cfg$umap_dims[cfg$umap_dims <= ncol(Embeddings(obj, "pca"))]
  message("Running UMAP from existing PCA dimensions: ", paste(range(dims), collapse = "-"))
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = dims, reduction.name = "umap", verbose = FALSE)
} else {
  warning("No PCA reduction found; skipping UMAP regeneration.")
}

message("Saving cleaned object: ", cfg$clean_major_rds)
make_dir(dirname(cfg$clean_major_rds))
saveRDS(obj, cfg$clean_major_rds)

DefaultAssay(obj) <- pick_assay(obj)
if (DefaultAssay(obj) == "SCT") {
  obj <- PrepSCTFindMarkers(obj, verbose = FALSE)
}

message("Finding major-cluster markers.")
markers <- FindAllMarkers(obj, min.pct = 0.1, only.pos = TRUE, verbose = FALSE)
save_csv(markers, file.path(output_dir, "major_cluster_markers.csv"))

counts <- obj@meta.data |>
  count(.data[[cfg$major_cluster_column]], .data[[cfg$sample_column]], name = "n_cells")
save_csv(counts, file.path(output_dir, "major_cluster_counts_by_sample.csv"))

if ("umap" %in% names(obj@reductions)) {
  p <- DimPlot(obj, group.by = cfg$major_cluster_column, label = FALSE, raster = FALSE) +
    scale_color_manual(values = cfg$major_colors, drop = FALSE) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    )
  save_plot(p, file.path(figure_dir, "major_cluster_umap.png"), width = 8, height = 6)
}

write_session_info(file.path(output_dir, "session_info.txt"))
message("Done.")
