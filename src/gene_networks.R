#!/usr/bin/env Rscript

# gene_networks.R
# Build a network for each age group of each tissue.

base_dir    <- normalizePath("../src")
data_dir    <- file.path(base_dir, "../data/normalized")
mapping_dir <- file.path(base_dir, "../data/raw")

# Create the top-level "mapped" directory if it doesn't exist
mapped_base_dir <- file.path(base_dir, "../data/mapped")
if (!dir.exists(mapped_base_dir)) {
  dir.create(mapped_base_dir, recursive = TRUE)
}

# Configurable threshold value
threshold_val <- 0.3

# Mapping file is located in data/raw, tab-separated:
# - First column: Gene symbol A
# - Second column: Gene symbol B
mapping_file <- "trrust_rawdata.human.tsv"
mapping_path <- file.path(mapping_dir, mapping_file)

# Read the mapping file, ignoring any columns beyond the first two
mapping_data <- read.delim(mapping_path, 
                           header = FALSE, 
                           sep = "\t", 
                           stringsAsFactors = FALSE)

colnames(mapping_data)[1:2] <- c("GeneA", "GeneB")
mapping_data <- mapping_data[, c("GeneA", "GeneB")]

# Helper function to get the short folder suffix from an age group like "20-29" => "20"
age_group_to_suffix <- function(col_name) {
  # Example: "20-29" -> "20", "30-39" -> "30", etc.
  strsplit(col_name, "-")[[1]][1]
}

# Find all normalized_{tissue name}.csv files in data/normalized
norm_files <- list.files(data_dir, pattern = "^normalized_.*\\.csv$", full.names = TRUE)
cat("Found", length(norm_files), "normalized tissue files in", data_dir, "\n")

for (nf in seq_along(norm_files)) {
  norm_file_path <- norm_files[nf]
  norm_file_name <- basename(norm_file_path)
  
  # Extract tissue name from "normalized_{tissue}.csv"
  tissue_name <- sub("^normalized_(.*)\\.csv$", "\\1", norm_file_name)
  
  cat(sprintf("\n[%d/%d] Processing tissue: %s\n", 
              nf, length(norm_files), tissue_name))
  
  # Read the normalized file
  norm_data <- read.csv(norm_file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # The first column is "Gene Symbols", the rest are age groups like "20-29", "30-39", ...
  age_group_cols <- colnames(norm_data)[-1]  # columns 2..end
  gene_symbols_set <- norm_data[["Gene Symbols"]]
  
  # Loop over each age group column
  for (ag_col in age_group_cols) {
    suffix <- age_group_to_suffix(ag_col)  # e.g. "20", "30", etc.
    
    # Create the subfolder in data/mapped/{suffix} (e.g. data/mapped/20)
    mapped_folder <- file.path(mapped_base_dir, suffix)
    if (!dir.exists(mapped_folder)) {
      dir.create(mapped_folder, recursive = TRUE)
    }
    
    # Output filename: mapped_{tissue name}_{age group}.csv
    # e.g. mapped_brain_cortex_20.csv
    output_file <- file.path(mapped_folder, 
                             paste0("mapped_", tissue_name, "_", suffix, ".csv"))
    
    # We'll collect edges in a list
    edges_list <- list()
    edge_idx   <- 1
    
    # This vector has expression values for the current age group
    expr_vector <- norm_data[[ag_col]]
    names(expr_vector) <- norm_data[["Gene Symbols"]]
    
    # Iterate over the TRRUST mappings
    for (row_i in seq_len(nrow(mapping_data))) {
      geneA <- mapping_data$GeneA[row_i]
      geneB <- mapping_data$GeneB[row_i]
      
      # 1) Check if both genes exist in the normalized data
      if (!(geneA %in% gene_symbols_set && geneB %in% gene_symbols_set)) {
        next
      }
      # 2) Check if geneA's expression in this age group > threshold_val
      valA <- expr_vector[geneA]
      if (is.na(valA) || valA <= threshold_val) {
        next
      }
      # If it passes, we add it to edges_list
      edges_list[[edge_idx]] <- c(geneA, geneB)
      edge_idx <- edge_idx + 1
    }
    
    # Convert edges_list to a data frame
    if (length(edges_list) == 0) {
      mapped_df <- data.frame(GeneA = character(0), 
                              GeneB = character(0),
                              stringsAsFactors = FALSE)
    } else {
      edges_mat <- do.call(rbind, edges_list)
      mapped_df <- data.frame(GeneA = edges_mat[,1],
                              GeneB = edges_mat[,2],
                              stringsAsFactors = FALSE)
    }
    
    # Write to CSV
    write.csv(mapped_df,
              file = output_file,
              row.names = FALSE,
              quote = FALSE)
    
    cat(sprintf("  Age group '%s': wrote %d edges -> %s\n", 
                ag_col, nrow(mapped_df), output_file))
  }
}

cat("\nAll network files have been created in data/mapped!\n")
