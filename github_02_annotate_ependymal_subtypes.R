#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

source("github_utils.R")
cfg <- load_config()
set.seed(cfg$random_seed)

output_dir <- make_dir(file.path(cfg$result_dir, "02_ependymal_subtypes"))
figure_dir <- make_dir(file.path(output_dir, "figures"))

message("Loading ependymal object: ", cfg$input_ependymal_rds)
ep <- read_seurat(cfg$input_ependymal_rds)
require_columns(ep@meta.data, cfg$ependymal_cluster_column, "Ependymal metadata")

cluster_id <- as.character(ep@meta.data[[cfg$ependymal_cluster_column]])
ep$ependymal_subtype <- NA_character_
for (subtype in names(cfg$ependymal_subtype_map)) {
  ep$ependymal_subtype[cluster_id %in% cfg$ependymal_subtype_map[[subtype]]] <- subtype
}
ep$ependymal_subtype <- factor(ep$ependymal_subtype, levels = names(cfg$ependymal_subtype_map))
Idents(ep) <- "ependymal_subtype"

message("Saving annotated ependymal object: ", cfg$annotated_ependymal_rds)
make_dir(dirname(cfg$annotated_ependymal_rds))
saveRDS(ep, cfg$annotated_ependymal_rds)

counts <- ep@meta.data |>
  count(ependymal_subtype, .data[[cfg$ependymal_cluster_column]], name = "n_cells")
save_csv(counts, file.path(output_dir, "ependymal_subtype_counts.csv"))

DefaultAssay(ep) <- pick_assay(ep)
if (DefaultAssay(ep) == "SCT") {
  ep <- PrepSCTFindMarkers(ep, verbose = FALSE)
}

message("Finding ependymal subtype markers.")
markers <- FindAllMarkers(ep, min.pct = 0.1, verbose = FALSE)
save_csv(markers, file.path(output_dir, "ependymal_subtype_markers.csv"))

if ("umap" %in% names(ep@reductions)) {
  p_umap <- DimPlot(ep, group.by = "ependymal_subtype", cols = cfg$ependymal_colors) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    )
  save_plot(p_umap, file.path(figure_dir, "ependymal_subtype_umap.png"), width = 6, height = 4)
}

genes <- intersect(cfg$marker_genes, rownames(ep))
if (length(genes) > 1) {
  avg <- AverageExpression(ep, assays = DefaultAssay(ep), features = genes, return.seurat = TRUE)
  mat <- GetAssayData(avg, assay = DefaultAssay(ep), layer = "scale.data")[genes, , drop = FALSE]
  write.csv(mat, file.path(output_dir, "ependymal_marker_average_expression.csv"), quote = FALSE)

  pdf(file.path(figure_dir, "ependymal_marker_heatmap.pdf"), width = 3.5, height = 6)
  pheatmap(
    mat,
    color = colorRampPalette(c("#5C88DA", "white", "#CC0C00"))(101),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = NA
  )
  dev.off()
}

write_session_info(file.path(output_dir, "session_info.txt"))
message("Done.")
