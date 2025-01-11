#!/usr/bin/env Rscript

# normalize_all_tissues.R
#
# This script looks in ../data/raw for any file matching the pattern:
#   gene_tpm_{tissue_name}.gctannotated.txt
# and for each one:
#   1) Reads the data
#   2) Removes the top epsilon% of gene expression
#   3) Min-max normalizes
#   4) Averages by age group
#   5) Saves the result to ../data/normalized/normalized_{tissue_name}.csv
#
# It reuses the logic from normalize_tissue.R, but loops over all matching files.
# Because the files are large, this may take a while, so we also print progress.

# -- Configurable parameters --------------------------------------------------

epsilon         <- 0.05    # top 5% removal
remove_zeroes   <- TRUE    # if TRUE, remove genes that end up with zero range
range_epsilon   <- 1e-12   # tiny tolerance for deciding if range is effectively zero
all_age_groups  <- c("20-29","30-39","40-49","50-59","60-69","70-79")

# -- Directory setup ---------------------------------------------------------

base_dir   <- normalizePath("../src")
input_dir  <- file.path(base_dir, "../data/raw")
output_dir <- file.path(base_dir, "../data/normalized")

# Make sure the output directory exists (create if needed)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -- Helper function: normalize_one_tissue -----------------------------------
normalize_one_tissue <- function(input_file_path, 
                                 epsilon = 0.05, 
                                 remove_zeroes = TRUE, 
                                 range_epsilon = 1e-12,
                                 age_groups = c("20-29","30-39","40-49","50-59","60-69","70-79")) {
  # 1) Determine how many gene columns exist
  line_4     <- readLines(input_file_path, n = 4)[4]
  total_cols <- length(strsplit(line_4, "\t")[[1]])
  gene_count <- total_cols - 4
  
  # 2) Read line 3 to get the gene symbols
  line_3       <- readLines(input_file_path, n = 3)[3]
  gene_symbols <- strsplit(line_3, "\t")[[1]]
  gene_symbols <- gene_symbols[1:gene_count]  # ignore metadata columns
  
  # 3) Read the expression + metadata (starting from line 4)
  lines_data <- read.delim(
    input_file_path,
    skip    = 3,
    header  = FALSE,
    sep     = "\t",
    quote   = "",
    stringsAsFactors = FALSE
  )
  
  # The last 4 columns are metadata, with the second-to-last as age group
  colnames(lines_data) <- c(
    gene_symbols,
    "meta1", 
    "meta2", 
    "meta3",  # second-to-last => age group
    "meta4"
  )
  
  # 4) Loop over each gene, compute normalized means by age group
  res_list <- vector("list", length = gene_count)
  row_idx  <- 1
  
  for (i in seq_len(gene_count)) {
    
    gene_name   <- gene_symbols[i]
    gene_values <- lines_data[[gene_name]]
    age_groups_ <- lines_data[["meta3"]]  # the age-group column
    
    # Remove top epsilon% of items
    cutoff      <- quantile(gene_values, probs = 1 - epsilon, na.rm = TRUE)
    trimmed_idx <- which(gene_values <= cutoff)
    
    trimmed_vals <- gene_values[trimmed_idx]
    trimmed_ages <- age_groups_[trimmed_idx]
    
    # Min/Max
    min_val   <- min(trimmed_vals, na.rm = TRUE)
    max_val   <- max(trimmed_vals, na.rm = TRUE)
    range_val <- max_val - min_val
    
    if (range_val < range_epsilon) {
      # If effectively zero range
      if (remove_zeroes) {
        # Skip storing
        next
      } else {
        # Otherwise, all values become 0
        normalized_vals <- rep(0, length(trimmed_vals))
      }
    } else {
      normalized_vals <- (trimmed_vals - min_val) / range_val
    }
    
    # Calculate average for each expected age group
    gene_means <- numeric(length(age_groups))
    for (ag_idx in seq_along(age_groups)) {
      ag     <- age_groups[ag_idx]
      idx_ag <- which(trimmed_ages == ag)
      
      if (length(idx_ag) == 0) {
        gene_means[ag_idx] <- NA
      } else {
        gene_means[ag_idx] <- mean(normalized_vals[idx_ag], na.rm = TRUE)
      }
    }
    
    # Store row
    row_vector         <- c(gene_name, gene_means)
    res_list[[row_idx]] <- row_vector
    row_idx            <- row_idx + 1
  }
  
  # Remove any unused list slots
  res_list <- res_list[!sapply(res_list, is.null)]
  
  if (length(res_list) == 0) {
    # Edge case: if all genes removed, create an empty data frame
    final_df <- data.frame(
      "Gene Symbols" = character(0),
      matrix(ncol=length(age_groups), nrow=0),
      check.names = FALSE
    )
    colnames(final_df) <- c("Gene Symbols", age_groups)
  } else {
    # Combine into a matrix, then convert to data frame
    final_mat  <- do.call(rbind, res_list)
    final_df   <- as.data.frame(final_mat, stringsAsFactors = FALSE)
    colnames(final_df) <- c("Gene Symbols", age_groups)
  }
  
  return(final_df)
}

# -- Main script logic: Loop over all files that match the pattern -----------
all_files <- list.files(input_dir, pattern = "^gene_tpm_.*\\.gctannotated\\.txt$", full.names = TRUE)

cat("Found", length(all_files), "raw tissue files.\n")

for (f in seq_along(all_files)) {
  file_path   <- all_files[f]
  
  # Extract the tissue name from the file name, e.g. gene_tpm_BRAIN_CORTEX.gctannotated.txt
  tissue_name <- sub("^gene_tpm_(.*)\\.gctannotated\\.txt$", "\\1", basename(file_path))
  
  # Print progress
  cat(sprintf("\n[%d/%d] Normalizing tissue '%s'...\n", f, length(all_files), tissue_name))
  
  # Perform normalization
  result_df <- normalize_one_tissue(
    input_file_path = file_path,
    epsilon         = epsilon,
    remove_zeroes   = remove_zeroes,
    range_epsilon   = range_epsilon,
    age_groups      = all_age_groups
  )
  
  # Construct output filename
  output_file <- file.path(output_dir, paste0("normalized_", tissue_name, ".csv"))
  
  # Write to CSV
  write.csv(
    result_df,
    file      = output_file,
    row.names = FALSE,
    quote     = FALSE
  )
  
  cat(sprintf("'%s' has been normalized and saved to %s\n", tissue_name, output_file))
}

cat("\nAll tissues have been processed!\n")
