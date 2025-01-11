#!/usr/bin/env Rscript

# normalize_tissue.R
# 
# Reads a .txt file of gene expression data:
#  - Tabs-separated
#  - Line 3: gene symbols
#  - Line 4 onward: each row = gene expression for subjects + 4 metadata columns
#     * second-to-last of those metadata columns = Age group (e.g. 20-29, 30-39, ..., 70-79)
#
# Steps:
#  1) Remove top epsilon% of each gene's expression values.
#  2) Min-max normalize the remaining values.
#  3) If remove_zeroes=TRUE, skip any genes whose range is effectively 0 after trimming.
#  4) Aggregate (mean) by age group.
#  5) Write results to normalized_<tissue_name>.csv

# -- Set up directories ------------------------------------------------------
base_dir    <- normalizePath("../src")
input_dir   <- file.path(base_dir, "../data/raw")
output_dir  <- file.path(base_dir, "../data/normalized")

# -- Input and output file details -------------------------------------------
input_file  <- "gene_tpm_brain_cortex.gctannotated.txt"
tissue_name <- sub("gene_tpm_(.*)\\.gctannotated\\.txt", "\\1", basename(input_file))
output_file <- file.path(output_dir, paste0("normalized_", tissue_name, ".csv"))
input_file_path <- file.path(input_dir, input_file)

# -- Configurable parameters -------------------------------------------------
epsilon         <- 0.05         # Top 5% removal
remove_zeroes   <- TRUE         # If TRUE, remove genes with effectively zero range
range_epsilon   <- 1e-12        # Floating-point tolerance for "zero" range

# -- Determine how many gene symbols there are -------------------------------
line_4     <- readLines(input_file_path, n = 4)[4]
total_cols <- length(strsplit(line_4, "\t")[[1]])
gene_count <- total_cols - 4
cat("Gene Symbols:", gene_count, "\n")

# -- Get gene symbols from line 3 --------------------------------------------
line_3       <- readLines(input_file_path, n = 3)[3]
gene_symbols <- strsplit(line_3, "\t")[[1]]
gene_symbols <- gene_symbols[1:gene_count]  # Ignore metadata columns
cat("Gene Symbols Stored\n")

# -- Read the expression + metadata (line 4 onward) into a data frame --------
lines_data <- read.delim(
  input_file_path,
  skip    = 3,        # Skip lines 1-3, so line 4 is row #1
  header  = FALSE,    # No header in these lines
  sep     = "\t",
  quote   = "",
  stringsAsFactors = FALSE
)

# lines_data has rows = subjects, columns = gene_count + 4
# The last 4 columns are metadata, with the second-to-last (index gene_count+3) being the age group.
colnames(lines_data) <- c(
  gene_symbols, 
  "meta1",    # 1st metadata col
  "meta2",    # 2nd metadata col
  "meta3",    # 3rd metadata col (SECOND to last => age group)
  "meta4"     # 4th metadata col (LAST)
)

# -- Predefine the expected age groups ---------------------------------------
all_age_groups <- c("20-29","30-39","40-49","50-59","60-69","70-79")

# We'll collect results in a list, then combine them at the end.
res_list <- vector("list", length = gene_count)
row_idx  <- 1

# -- Main loop: for each gene, normalize, then average by age group ----------
for (i in seq_len(gene_count)) {
  
  # 1) Get gene expression
  gene_name   <- gene_symbols[i]
  gene_values <- lines_data[[gene_name]]
  
  # 2) Get age groups
  age_groups  <- lines_data[["meta3"]]
  
  # 3) Remove top epsilon% (5%) of items
  cutoff      <- quantile(gene_values, probs = 1 - epsilon, na.rm = TRUE)
  trimmed_idx <- which(gene_values <= cutoff)
  
  trimmed_vals <- gene_values[trimmed_idx]
  trimmed_ages <- age_groups[trimmed_idx]
  
  # 4) Min-max normalization (with range epsilon check)
  min_val   <- min(trimmed_vals, na.rm = TRUE)
  max_val   <- max(trimmed_vals, na.rm = TRUE)
  range_val <- max_val - min_val
  
  if (range_val < range_epsilon) {
    # If the range is effectively zero
    if (remove_zeroes) {
      # Skip storing this gene entirely
      next
    } else {
      # If not removing, all normalized values become 0
      normalized_vals <- rep(0, length(trimmed_vals))
    }
  } else {
    # Normal min-max
    normalized_vals <- (trimmed_vals - min_val) / range_val
  }
  
  # 5) Compute average for each age group
  gene_means <- numeric(length(all_age_groups))
  for (ag_idx in seq_along(all_age_groups)) {
    ag  <- all_age_groups[ag_idx]
    idx <- which(trimmed_ages == ag)
    
    if (length(idx) == 0) {
      gene_means[ag_idx] <- NA
    } else {
      gene_means[ag_idx] <- mean(normalized_vals[idx], na.rm = TRUE)
    }
  }
  
  # 6) Store the row in res_list
  row_vector      <- c(gene_name, gene_means)
  res_list[[row_idx]] <- row_vector
  row_idx         <- row_idx + 1
}

# -- Combine all stored rows into a final data frame -------------------------
res_list <- res_list[ !sapply(res_list, is.null) ]  # remove empty entries
if (length(res_list) == 0) {
  # Edge case: if all genes were removed, create an empty data frame
  final_df <- data.frame(
    "Gene Symbols" = character(0),
    matrix(ncol=length(all_age_groups), nrow=0),
    check.names = FALSE
  )
  colnames(final_df) <- c("Gene Symbols", all_age_groups)
} else {
  final_mat  <- do.call(rbind, res_list)
  final_df   <- as.data.frame(final_mat, stringsAsFactors = FALSE)
  colnames(final_df) <- c("Gene Symbols", all_age_groups)
}

# -- Write the result to CSV -------------------------------------------------
write.csv(
  final_df,
  file      = output_file,
  row.names = FALSE,
  quote     = FALSE
)

cat("Normalization complete!\n")
if (remove_zeroes) {
  cat("Genes with effectively zero range (<", range_epsilon, ") were removed.\n")
}
cat("Output written to:", output_file, "\n")
