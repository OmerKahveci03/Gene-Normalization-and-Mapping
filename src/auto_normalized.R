# This script takes raw data and converts it into the normalized version that we will be using # nolint
# Due to the size of the data, this will take several hours to process 49 different tissues # nolint
# Load necessary libraries
library(data.table)
library(dplyr)
library(openxlsx)
library(tidyr)

# Define file paths
base_dir <- normalizePath("../")
input_dir <- file.path(base_dir, "data/raw")
output_dir <- file.path(base_dir, "data/normalized")

# Get a list of all gene_tpm files in the input directory
input_files <- list.files(input_dir, pattern = "gene_tpm_.*\\.gctannotated\\.txt", full.names = TRUE) # nolint

# Define a function to process each file
process_file <- function(file_path) {
  # Extract tissue name from the file name
  tissue_name <- sub("gene_tpm_(.*)\\.gctannotated\\.txt", "\\1", basename(file_path)) # nolint
  # Define output file paths
  output_file_name <- paste0("normalized_", tissue_name, ".txt")
  output_file_path <- file.path(output_dir, output_file_name)
  # Get expression count
  line_4 <- readLines(file_path, n = 4)[4]
  expression_counts <- length(strsplit(line_4, "\t")[[1]]) - 4
  cat("Gene Symbols:", expression_counts, "\n")
  # Get gene symbols
  line_3 <- readLines(file_path, n = 3)[3]
  gene_symbols <- strsplit(line_3, "\t")[[1]]
  gene_symbols <- gene_symbols[1:expression_counts]
  cat("Gene Symbols Stored", "\n")
  # Filter out gene symbols with underscores
  valid_gene_symbols <- gene_symbols[!grepl("_", gene_symbols)]
  expression_counts <- length(valid_gene_symbols)
  # Ensure unique gene symbols
  valid_gene_symbols <- make.unique(valid_gene_symbols)
  # Initialize Data Storage
  data_lines <- list()
  # Read the expression data starting from the 4th line
  con <- file(file_path, "r")
  open(con)
  for (i in 1:3) { readLines(con, n = 1) }  # Skip the first 3 lines
  while (length(line <- readLines(con, n = 1)) > 0) {
    if (nchar(line) > 0) {
      # Split the line by tab and keep only the required number of elements
      line_data <- unlist(strsplit(line, "\t"))
      # Extract the age group
      age_group <- line_data[(length(line_data) - 1)]
      # Extract the gene expression data
      expression_data <- line_data[1:expression_counts]
      # Combine age group and expression data
      combined_data <- c(Age = age_group, expression_data)
      # Append the extracted data to the list
      data_lines <- c(data_lines, list(combined_data))
    }
  }
  close(con)
  # Convert the list of data lines to a data frame
  data_df <- as.data.frame(do.call(rbind, data_lines), stringsAsFactors = FALSE)
  colnames(data_df) <- c("Age", valid_gene_symbols)
  # Convert expression columns to numeric
  data_df[valid_gene_symbols] <- lapply(data_df[valid_gene_symbols], as.numeric)
  # Filter out gene symbols where max - min < epsilon before normalization
  epsilon <- 0.01
  max_min_diff <- data_df %>%
    select(-Age) %>%
    summarise(across(everything(), ~ max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE))) #nolint
  valid_gene_symbols <- names(max_min_diff)[max_min_diff > epsilon]
  filtered_data_df <- data_df %>% select(Age, all_of(valid_gene_symbols))
  # Normalize the filtered gene expression data
  normalized_data <- filtered_data_df %>%
    select(-Age) %>%
    mutate(across(everything(), ~ (.) / (max(., na.rm = TRUE) - min(., na.rm = TRUE)))) #nolint
  # Add Age column back to the normalized data
  normalized_data <- cbind(Age = filtered_data_df$Age, normalized_data)
  # Normalize gene symbols by removing parts after the dash
  normalized_data <- normalized_data %>%
    pivot_longer(cols = -Age, names_to = "GeneSymbol", values_to = "Expression") %>% #nolint
    mutate(GeneSymbol = sub("-.*", "", GeneSymbol))
  # Group by Age and normalized Gene Symbols, and calculate the average expression for each group #nolint
  avg_expression_df <- normalized_data %>%
    group_by(Age, GeneSymbol) %>%
    summarise(AvgExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
  # Reshape the data to wide format
  avg_expression_df_wide <- avg_expression_df %>%
    pivot_wider(names_from = GeneSymbol, values_from = AvgExpression)
  # Sort the data frame by Age groups
  age_order <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99") #nolint
  avg_expression_df_wide$Age <- factor(avg_expression_df_wide$Age, levels = age_order) #nolint
  avg_expression_df_wide <- avg_expression_df_wide %>% arrange(Age) %>% na.omit() #nolint
  # Write the final product to a text file
  write.table(avg_expression_df_wide, output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE) #nolint
  # Write the final product to an Excel file
  write.xlsx(avg_expression_df_wide, output_excel_path)
  cat("Processed:", file_path, "\n")
}

# Process each file in the input directory
for (file_path in input_files) {
  process_file(file_path)
}

cat("Reached End", "\n")
