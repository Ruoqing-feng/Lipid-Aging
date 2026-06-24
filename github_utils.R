#!/usr/bin/env Rscript

load_config <- function(default_file = "github_config.R") {
  args <- commandArgs(trailingOnly = TRUE)
  config_arg <- grep("^--config=", args, value = TRUE)
  config_file <- if (length(config_arg) == 1) {
    sub("^--config=", "", config_arg)
  } else {
    default_file
  }

  if (!file.exists(config_file)) {
    stop(
      "Config file not found: ", config_file,
      "\nCopy github_config.example.R to github_config.R and edit it.",
      call. = FALSE
    )
  }

  env <- new.env(parent = globalenv())
  sys.source(config_file, envir = env)
  env
}

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

require_columns <- function(data, columns, label = "data") {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0) {
    stop(label, " is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

read_seurat <- function(path) {
  if (!file.exists(path)) {
    stop("Input RDS not found: ", path, call. = FALSE)
  }
  readRDS(path)
}

pick_assay <- function(object, preferred = c("SCT", "RNA")) {
  assays <- names(object@assays)
  hit <- preferred[preferred %in% assays]
  if (length(hit) > 0) {
    hit[1]
  } else {
    DefaultAssay(object)
  }
}

save_csv <- function(data, path) {
  make_dir(dirname(path))
  write.csv(data, path, row.names = FALSE, quote = FALSE)
  invisible(path)
}

save_plot <- function(plot, path, width = 7, height = 5, dpi = 300) {
  make_dir(dirname(path))
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
  invisible(path)
}

repair_legacy_spatial_slots <- function(object) {
  for (image_name in names(object@images)) {
    image_object <- object@images[[image_name]]
    if (!inherits(image_object, "FOV")) {
      next
    }

    if (inherits(try(methods::slot(image_object, "misc"), silent = TRUE), "try-error")) {
      methods::slot(image_object, "misc") <- list()
    }

    for (boundary_name in names(image_object@boundaries)) {
      boundary_object <- image_object@boundaries[[boundary_name]]
      if (
        inherits(boundary_object, "Segmentation") &&
          inherits(try(methods::slot(boundary_object, "sf.data"), silent = TRUE), "try-error")
      ) {
        methods::slot(boundary_object, "sf.data") <- NULL
        image_object@boundaries[[boundary_name]] <- boundary_object
      }
    }

    object@images[[image_name]] <- image_object
  }
  object
}

write_session_info <- function(path) {
  make_dir(dirname(path))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}
