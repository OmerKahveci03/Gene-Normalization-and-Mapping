#!/usr/bin/env Rscript

# age_based_influence_subset.R
# Only generates plots for: 
#   • max → Thyroid 
#   • min → Prostate

# Libraries & Utilities
library(igraph)
library(readxl)
source("utility.R")  # must provide compute_influence_list() and age_based_influence()

# Configurable values
lambda     <- 1
epsilon    <- 0.005
age_groups <- c(20, 30, 40, 50, 60, 70)

# Paths
project_root <- normalizePath("..")
excel_dir    <- file.path(project_root, "data", "age_indexes", "excel_files")
max_file     <- file.path(excel_dir, "a1_max_age_indexes.xlsx")
min_file     <- file.path(excel_dir, "a1_min_age_indexes.xlsx")

# Output directories (will be created if missing)
out_max_dir <- file.path(project_root, "visuals/final")
out_min_dir <- file.path(project_root, "visuals/final")
if (!dir.exists(out_max_dir)) dir.create(out_max_dir, recursive = TRUE)
if (!dir.exists(out_min_dir)) dir.create(out_min_dir, recursive = TRUE)

# Helper: read 'best' sheet and parse Gene_Set into a list
read_best_sheet <- function(path) {
  df <- read_excel(path, sheet = "best")
  # Expect columns: Tissue | Metric | Gene_Set
  df$Gene_Set <- strsplit(df$Gene_Set, ",\\s*")
  df
}

# Plotting function (identical to your original)
plot_for_tissue <- function(tissue, genes, type) {
  start_time <- Sys.time()
  cat(sprintf("[%s] Starting %s for tissue '%s'\n",
              format(start_time, "%Y-%m-%d %H:%M:%S"),
              toupper(type), tissue))
  
  # Compute influence
  influence_list <- compute_influence_list(
    age_groups  = age_groups,
    tissue_name = tissue,
    gene_set    = genes,
    lambda      = lambda,
    data_dir    = project_root
  )
  
  # Log per-age values
  for (k in seq_along(age_groups)) {
    cat(sprintf("  Age %2d => Influence: %f\n",
                age_groups[k], influence_list[k]))
  }
  # Overall age-based influence
  abi_val <- age_based_influence(influence_list, epsilon)
  cat(sprintf("  Age-based influence: %f\n", abi_val))
  
  # Y-axis limit
  max_val <- max(influence_list, na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 0.1
  ymax <- ceiling(max_val * 10) / 10
  
  # Plot
  out_dir <- if (type == "max") out_max_dir else out_min_dir
  file_name <- sprintf("%s_%s_age_index.png", type, tolower(tissue))
  png(filename = file.path(out_dir, file_name), width = 800, height = 600)
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(
    x    = age_groups,
    y    = influence_list,
    type = "l",
    pch  = 19,
    col  = if (type == "max") "darkgreen" else "darkred",
    lwd = 2.3,
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
    xlab = "Age Group",
    ylab = "Influence(S, V)",
    main = sprintf("%s Age Index: %s", toupper(type), tissue),
    ylim = c(0, ymax)
  )
  grid()
  dev.off()
  
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf("[%s] Finished %s for tissue '%s' (%.2f sec)\n\n",
              format(end_time, "%Y-%m-%d %H:%M:%S"),
              toupper(type), tissue, duration))
}


# --- MAIN ---

# Read only the 'best' sheets
max_df <- read_best_sheet(max_file)
min_df <- read_best_sheet(min_file)

# Filter down to Thyroid in max, Prostate in min
max_thyroid <- subset(max_df, tolower(Tissue) == "thyroid")
min_prostate <- subset(min_df, tolower(Tissue) == "prostate")

# Plot Thyroid (max)
if (nrow(max_thyroid) == 1) {
  plot_for_tissue(
    tissue = max_thyroid$Tissue,
    genes  = max_thyroid$Gene_Set[[1]],
    type   = "max"
  )
} else {
  stop("Thyroid not found or multiple entries in max sheet")
}

# Plot Prostate (min)
if (nrow(min_prostate) == 1) {
  plot_for_tissue(
    tissue = min_prostate$Tissue,
    genes  = min_prostate$Gene_Set[[1]],
    type   = "min"
  )
} else {
  stop("Prostate not found or multiple entries in min sheet")
}

cat("Done. Plots are in:\n", out_max_dir, "\n", out_min_dir, "\n")
