#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("github_utils.R")
cfg <- load_config()

output_dir <- make_dir(file.path(cfg$result_dir, "03_spatial_neighborhood_summary"))
figure_dir <- make_dir(file.path(output_dir, "figures"))

message("Loading cleaned major object: ", cfg$clean_major_rds)
major <- read_seurat(cfg$clean_major_rds)
message("Loading annotated ependymal object: ", cfg$annotated_ependymal_rds)
ep <- read_seurat(cfg$annotated_ependymal_rds)

required_meta <- c(
  cfg$major_cluster_column, cfg$sample_column, cfg$age_column,
  cfg$position_column, cfg$x_column, cfg$y_column
)
require_columns(major@meta.data, required_meta, "Major-object metadata")
require_columns(ep@meta.data, "ependymal_subtype", "Ependymal metadata")

meta <- major@meta.data
meta$cell <- rownames(meta)
meta$ependymal_subtype <- NA_character_
matched <- match(Cells(ep), meta$cell)
meta$ependymal_subtype[matched[!is.na(matched)]] <- as.character(ep$ependymal_subtype)[!is.na(matched)]

sample_col <- cfg$sample_column
x_col <- cfg$x_column
y_col <- cfg$y_column
cluster_col <- cfg$major_cluster_column
radius <- cfg$neighbor_radius

message("Calculating neighbors within radius: ", radius)
neighbor_rows <- vector("list", sum(!is.na(meta$ependymal_subtype)))
row_i <- 1

for (sample_id in sort(unique(meta[[sample_col]]))) {
  sample_meta <- meta[meta[[sample_col]] == sample_id, , drop = FALSE]
  ep_meta <- sample_meta[!is.na(sample_meta$ependymal_subtype), , drop = FALSE]
  if (nrow(ep_meta) == 0) {
    next
  }

  coords <- as.matrix(sample_meta[, c(x_col, y_col)])
  for (i in seq_len(nrow(ep_meta))) {
    center <- as.numeric(ep_meta[i, c(x_col, y_col)])
    distance <- sqrt((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2)
    is_neighbor <- distance > 0 & distance <= radius
    neighbor_table <- table(sample_meta[[cluster_col]][is_neighbor])

    neighbor_rows[[row_i]] <- data.frame(
      ependymal_cell = ep_meta$cell[i],
      sample = sample_id,
      age = ep_meta[[cfg$age_column]][i],
      position = ep_meta[[cfg$position_column]][i],
      ependymal_subtype = ep_meta$ependymal_subtype[i],
      neighbor_major_cluster = names(neighbor_table),
      n_neighbors = as.integer(neighbor_table),
      radius = radius
    )
    row_i <- row_i + 1
  }
}

neighbors <- bind_rows(neighbor_rows)
save_csv(neighbors, file.path(output_dir, "ependymal_neighbor_major_clusters.csv"))

summary <- neighbors |>
  group_by(age, position, ependymal_subtype, neighbor_major_cluster) |>
  summarize(mean_neighbors = mean(n_neighbors), .groups = "drop")
save_csv(summary, file.path(output_dir, "neighbor_summary_by_age_position.csv"))

if (nrow(summary) > 0) {
  p <- ggplot(summary, aes(x = neighbor_major_cluster, y = mean_neighbors, fill = ependymal_subtype)) +
    geom_col(position = "dodge") +
    facet_grid(age ~ position, scales = "free_y") +
    scale_fill_manual(values = cfg$ependymal_colors, drop = FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    ) +
    labs(x = NULL, y = "Mean neighboring cells")
  save_plot(p, file.path(figure_dir, "neighbor_summary_by_age_position.png"), width = 10, height = 6)
}

write_session_info(file.path(output_dir, "session_info.txt"))
message("Done.")
