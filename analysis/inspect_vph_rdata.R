#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: inspect_vph_rdata.R <input.rdata> <output-dir>")
}
input <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

env <- new.env(parent = baseenv())
loaded <- load(input, envir = env)

rows <- list()
for (name in loaded) {
  obj <- get(name, envir = env)
  d <- dim(obj)
  rows[[length(rows) + 1]] <- data.frame(
    object = name,
    class = paste(class(obj), collapse = ";"),
    nrow = if (is.null(d)) NA_integer_ else d[[1]],
    ncol = if (is.null(d) || length(d) < 2) NA_integer_ else d[[2]],
    size_bytes = as.numeric(object.size(obj)),
    stringsAsFactors = FALSE
  )
}
inv <- do.call(rbind, rows)
write.csv(inv, file.path(outdir, "object_inventory.csv"), row.names = FALSE)

preview_lines <- c("# VPH RData object inspection", "")
for (name in loaded) {
  obj <- get(name, envir = env)
  preview_lines <- c(preview_lines, paste0("## `", name, "`"))
  preview_lines <- c(preview_lines, paste0("- class: `", paste(class(obj), collapse = ";"), "`"))
  d <- dim(obj)
  if (!is.null(d)) preview_lines <- c(preview_lines, paste0("- dim: ", paste(d, collapse = " x ")))
  if (!is.null(colnames(obj))) {
    preview_lines <- c(preview_lines, paste0("- columns/head: `", paste(head(colnames(obj), 30), collapse = "`, `"), "`"))
  }
  if (!is.null(rownames(obj))) {
    preview_lines <- c(preview_lines, paste0("- rownames/head: `", paste(head(rownames(obj), 10), collapse = "`, `"), "`"))
  }
  preview_lines <- c(preview_lines, "")
}
writeLines(preview_lines, file.path(outdir, "object_schema.md"))

# Candidate metadata tables: modest column count and cell-like row count.
meta_candidates <- c()
for (name in loaded) {
  obj <- get(name, envir = env)
  if ((is.data.frame(obj) || is.matrix(obj)) && !is.null(dim(obj))) {
    nr <- nrow(obj); nc <- ncol(obj)
    cols <- tolower(colnames(obj) %||% character())
    keyword <- any(grepl("cluster|class|sex|sample|batch|umap|chem", cols))
    cellshape <- nr >= 1000 && nr <= 100000 && nc <= 100
    if (keyword || cellshape) {
      meta_candidates <- c(meta_candidates, name)
      if (nc <= 100) {
        out <- as.data.frame(obj)
        out$.rowname <- rownames(obj)
        write.csv(out, file.path(outdir, paste0("candidate_", gsub("[^A-Za-z0-9_.-]", "_", name), ".csv")), row.names = FALSE)
      }
    }
  }
}
writeLines(meta_candidates, file.path(outdir, "metadata_candidates.txt"))

cat("Loaded objects:\n")
print(inv)
cat("\nMetadata candidates:\n")
print(meta_candidates)
