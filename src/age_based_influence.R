# age_based_influence.R

# Libraries & Utilities
library(igraph)
library(readxl)
source("utility.R")  # Brings in compute_influence_list() and age_based_influence()

# Configurable Values
lambda       <- 1
epsilon      <- 0.005
age_groups   <- c(20, 30, 40, 50, 60, 70)

# Paths
project_root <- normalizePath("..")
excel_dir    <- file.path(project_root, "data", "age_indexes", "excel_files")
max_file     <- file.path(excel_dir, "a1_max_age_indexes.xlsx")
min_file     <- file.path(excel_dir, "a1_min_age_indexes.xlsx")

# Output directories
out_max_dir <- file.path(project_root, "visuals", "max_age_plots")
out_min_dir <- file.path(project_root, "visuals", "min_age_plots")
if (!dir.exists(out_max_dir)) dir.create(out_max_dir, recursive = TRUE)
if (!dir.exists(out_min_dir)) dir.create(out_min_dir, recursive = TRUE)

# Helper: read best sheet and parse gene set
read_best_sheet <- function(path) {
  df <- read_excel(path, sheet = "best")
  # Columns: Tissue, Metric, Gene_Set
  df$Gene_Set <- strsplit(df$Gene_Set, ",\\s*")
  return(df)
}

# Load both max and min tables
max_df <- read_best_sheet(max_file)
min_df <- read_best_sheet(min_file)

# Function to generate and save plots for a given tissue
plot_for_tissue <- function(tissue, genes, type) {
  start_time <- Sys.time()
  cat(sprintf("[%s] Starting %s for tissue '%s' at %s\n",
              format(start_time, "%Y-%m-%d %H:%M:%S"), toupper(type), tissue,
              format(start_time, "%H:%M:%S")))
  
  # Compute influence list across ages
  influence_list <- compute_influence_list(
    age_groups  = age_groups,
    tissue_name = tissue,
    gene_set    = genes,
    lambda      = lambda,
    data_dir    = project_root
  )
  
  # Log individual age-group values
  for (k in seq_along(age_groups)) {
    cat(sprintf("  Age %2d => Influence: %f\n", age_groups[k], influence_list[k]))
  }
  
  # Compute overall age-based influence value (optional)
  abi_val <- age_based_influence(influence_list, epsilon)
  cat(sprintf("  Age-based influence: %f\n", abi_val))
  
  # Determine y-axis limit
  max_val <- max(influence_list, na.rm = TRUE)
  if (is.infinite(max_val) || max_val <= 0) max_val <- 0.1
  ymax <- ceiling(max_val * 10) / 10
  
  # Create plot
  out_dir <- if (type == "max") out_max_dir else out_min_dir
  file_name <- sprintf("%s_%s_age_index.png", type, tissue)
  png(filename = file.path(out_dir, file_name), width = 800, height = 600)
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(
    x    = age_groups,
    y    = influence_list,
    type = "b",
    pch  = 19,
    col  = if (type == "max") "darkgreen" else "darkred",
    xlab = "Age Group",
    ylab = "Influence(S, V)",
    main = sprintf("%s Age Index: %s", toupper(type), tissue),
    ylim = c(0, ymax)
  )
  grid()
  dev.off()
  
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf("[%s] Finished %s for tissue '%s' at %s (%.2f sec)\n",
              format(end_time, "%Y-%m-%d %H:%M:%S"), toupper(type), tissue,
              format(end_time, "%H:%M:%S"), duration))
}

# Iterate over max tissues
for (i in seq_len(nrow(max_df))) {
  row <- max_df[i, ]
  plot_for_tissue(row$Tissue, row$Gene_Set[[1]], "max")
}

# Iterate over min tissues
for (i in seq_len(nrow(min_df))) {
  row <- min_df[i, ]
  plot_for_tissue(row$Tissue, row$Gene_Set[[1]], "min")
}

cat("All plots created in:", out_max_dir, "and", out_min_dir, "\n")
