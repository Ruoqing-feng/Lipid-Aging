# Ependymal Aging Spatial Transcriptomics

This repository contains reproducible R scripts for analyzing ependymal-cell aging patterns in spatial transcriptomics data. The workflow starts from Seurat objects, removes low-quality major-cell annotations, annotates ependymal subtypes, and summarizes local cellular neighborhoods around ependymal cells.

## Repository Layout

```text
Analysis/
  scripts/
    github_config.example.R
    github_utils.R
    github_01_clean_major_object.R
    github_02_annotate_ependymal_subtypes.R
    github_03_spatial_neighborhood_summary.R
  data/
    input Seurat RDS files, not tracked if large
  results/
    generated tables, plots, and session logs
```

## Requirements

Install R and the following packages:

```r
install.packages(c("dplyr", "ggplot2", "pheatmap"))
install.packages("Seurat")
```

The scripts expect Seurat objects with metadata columns for major cell type, sample, age, position, and spatial coordinates. Column names are configurable in `github_config.R`.

## Quick Start

From `Analysis/scripts`:

```bash
cp github_config.example.R github_config.R
Rscript github_01_clean_major_object.R --config=github_config.R
Rscript github_02_annotate_ependymal_subtypes.R --config=github_config.R
Rscript github_03_spatial_neighborhood_summary.R --config=github_config.R
```

Before running, edit `github_config.R` so the input paths point to your local Seurat `.rds` files.

## Workflow

1. `github_01_clean_major_object.R`
   - Loads a clustered Seurat object.
   - Removes cells assigned to the configured bad-cell label.
   - Regenerates UMAP from an existing PCA reduction when available.
   - Exports marker genes, sample-level cell counts, a UMAP plot, and a cleaned RDS object.

2. `github_02_annotate_ependymal_subtypes.R`
   - Loads an ependymal-cell Seurat object.
   - Maps clustering IDs to interpretable ependymal subtype labels.
   - Exports subtype counts, marker genes, UMAP, marker heatmap, and an annotated RDS object.

3. `github_03_spatial_neighborhood_summary.R`
   - Combines the cleaned major-cell object with annotated ependymal subtypes.
   - Counts neighboring major-cell classes around each ependymal cell within a configurable spatial radius.
   - Exports per-cell neighbor tables, grouped summaries, and a summary figure.

## Data

Large input and output objects are intentionally excluded from source control. Recommended tracked files are:

- R scripts under `Analysis/scripts/`
- Small example metadata tables, if available
- Documentation and configuration examples

Recommended ignored files are:

- `*.rds`
- large raw count matrices
- generated figures and result tables

## Reproducibility

Each script writes `session_info.txt` into its result folder. Include those logs when reporting results so package versions and runtime details can be checked later.

## Citation

If you use this workflow, please cite the related study or dataset once available.
