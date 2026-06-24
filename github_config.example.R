# Copy this file to github_config.R and edit paths before running the scripts.

project_dir <- normalizePath("..", mustWork = FALSE)
data_dir <- file.path(project_dir, "data")
result_dir <- file.path(project_dir, "results")

input_major_rds <- file.path(data_dir, "lipidall_clustered.rds")
clean_major_rds <- file.path(data_dir, "lipidall_clustered_major_clean.rds")
input_ependymal_rds <- file.path(data_dir, "lipidependy_clean.rds")
annotated_ependymal_rds <- file.path(data_dir, "lipidependy_clean_annotated.rds")

major_cluster_column <- "majorcluster"
bad_cluster_label <- "Badcell"
sample_column <- "sample"
age_column <- "Age"
position_column <- "position"
x_column <- "x"
y_column <- "y"

ependymal_cluster_column <- "SCT_snn_res.0.3"
ependymal_subtype_map <- list(
  ependymal1 = c("0", "1"),
  ependymal2 = c("4", "5")
)

major_colors <- c(
  Astro = "#E64B35",
  ChP = "#00A087",
  Endo = "#3C5488",
  Ependymal = "#F39B7F",
  Micro = "#8491B4",
  Neuron = "#91D1C2",
  Oligo = "#DC0000",
  Opc = "#7E6148"
)

ependymal_colors <- c(
  ependymal1 = "#D97706",
  ependymal2 = "#7E22CE"
)

marker_genes <- c(
  "Map2", "Rxrg", "Gpm6b", "Itgb1", "Kcnk2", "Plin4", "Rarres2",
  "Ttyh2", "Ldlr", "Chat", "Ifit1", "C4b", "Rsad2", "Isg15",
  "Abca1", "Abcg1", "Bst2", "H2-K1", "Ifitm3", "Gsn", "Dusp1"
)

umap_dims <- 1:40
neighbor_radius <- 100
random_seed <- 2026
