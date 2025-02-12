# drug_efficacy.R - takes roughly 9 minutes
# -----------------------------------------------------------
# Input:
#   - Disease genes S (from filtered_neurodegenerative_diseases.csv) for Alzheimer Disease.
#   - Drug candidate information from filtered_drug_gene_dataset.xlsx:
#         * "Sorted" sheet: Column A = drug name, Column B = gene count (sorted descending).
#         * "Raw" sheet: Column A = drug name, Column B = gene name (may be duplicated).
#
#   For each candidate drug (with at least 20 target genes), we use its full target gene set
#   to calculate efficacy across tissues.
#
#   The efficacy formula for a set of drug targets D is:
#      efficacy(D, S) = ( sum_{s in S'} influence_S_g(graph, D, s, lambda) ) /
#                       ( sum_{n in N'} influence_S_g(graph, D, n, lambda) )
#      where S' is the set of disease genes present in the network and N' is all other genes.
#
#   Output:
#     - A CSV file (data/drug/efficacy.csv) with tissues in rows and a column for each candidate drug.
#
library(igraph)
library(readxl)
source("utility.R")

nd_file <- file.path("..", "data", "drug", "filtered_neurodegenerative_diseases.csv")
if (!file.exists(nd_file)) {
  stop(sprintf("Neurodegenerative diseases file not found: %s", nd_file))
}

# Read the CSV file (adjust header parameter if needed)
nd_data <- read.csv(nd_file, header = FALSE, stringsAsFactors = FALSE)
# Assume Column 1 = gene name, Column 2 = disease name.
disease_genes <- nd_data[nd_data[, 2] == "Alzheimer Disease", 1]
if (length(disease_genes) == 0) {
  stop("No genes found for 'Alzheimer Disease' in the neurodegenerative diseases file.")
}
# Clean gene names.
disease_genes <- trimws(disease_genes)

drug_excel_file <- file.path("..", "data", "drug", "filtered_drug_gene_dataset.xlsx")
if (!file.exists(drug_excel_file)) {
  stop(sprintf("Drug Excel file not found: %s", drug_excel_file))
}

# Read the "Sorted" sheet: Column A = drug name, Column B = gene count.
sorted_data <- read_excel(drug_excel_file, sheet = "Sorted", col_names = FALSE)
# Keep only those drugs with a gene count >= 20.
candidate_sorted <- sorted_data[sorted_data[[2]] >= 20, ]
if(nrow(candidate_sorted) == 0) {
  stop("No candidate drugs with gene count >= 20 found in the Sorted sheet.")
}
# Extract candidate drug names.
candidate_drugs <- candidate_sorted[[1]]

# Read the "Raw" sheet: Column A = drug name, Column B = gene name.
raw_data <- read_excel(drug_excel_file, sheet = "Raw", col_names = FALSE)

# Build a named list of candidate drug targets.
# For each candidate drug, retrieve its target genes from the "Raw" sheet.
candidate_drug_targets_list <- lapply(candidate_drugs, function(drug) {
  targets <- raw_data[[2]][raw_data[[1]] == drug]
  unique(trimws(targets))
})
names(candidate_drug_targets_list) <- candidate_drugs

# Compute the intersection count for each candidate drug.
candidate_intersections <- sapply(candidate_drugs, function(drug) {
  targets <- candidate_drug_targets_list[[drug]]
  length(intersect(targets, disease_genes))
})
# Order candidate drugs in descending order of intersection count.
ordered_indices <- order(candidate_intersections, decreasing = TRUE)
# Pick the top n candidate drugs.
top_n <- 270
if(length(ordered_indices) < top_n) {
  top_n <- length(ordered_indices)
}
top_candidate_indices <- ordered_indices[1:top_n]
candidate_drugs <- candidate_drugs[top_candidate_indices]
candidate_drug_targets_list <- candidate_drug_targets_list[candidate_drugs]

cat("Top candidate drugs (with intersection counts):\n")
for(i in seq_along(candidate_drugs)) {
  cat(sprintf("  %s: %d\n", candidate_drugs[i], candidate_intersections[top_candidate_indices[i]]))
}

age_group <- 20
lambda    <- 1

mapped_dir <- file.path("..", "data", "mapped", as.character(age_group))
if (!dir.exists(mapped_dir)) {
  stop(sprintf("Mapped directory does not exist: %s", mapped_dir))
}

mapped_files <- list.files(
  path       = mapped_dir,
  pattern    = paste0("^mapped_.*_", age_group, "\\.csv$"),
  full.names = TRUE
)
if (length(mapped_files) == 0) {
  stop(sprintf("No mapped CSV files found in %s", mapped_dir))
}

drug_dir <- file.path("..", "data", "drug")
if (!dir.exists(drug_dir)) {
  dir.create(drug_dir, recursive = TRUE)
}
output_file <- file.path(drug_dir, "efficacy.csv")

calc_efficacy <- function(graph, drug_targets, disease_genes, non_disease_genes, lambda = 1) {
  # Only consider drug targets that are in the graph.
  valid_drug_targets <- intersect(drug_targets, V(graph)$name)
  if (length(valid_drug_targets) == 0) {
    return(0)
  }
  
  # Only consider disease genes that are in the graph.
  valid_disease_genes <- intersect(disease_genes, V(graph)$name)
  if (length(valid_disease_genes) == 0) {
    return(0)
  }
  
  # Precompute distances from the set of drug targets to all nodes in the graph.
  dist_matrix <- distances(graph, v = valid_drug_targets)
  
  # If multiple sources, take the minimum distance for each node.
  if (length(valid_drug_targets) > 1) {
    min_distances <- apply(dist_matrix, 2, min)
  } else {
    min_distances <- as.vector(dist_matrix)
  }
  
  # Compute influence values for all nodes.
  influence_values <- exp(-lambda * min_distances)
  
  # Calculate numerator: sum of influences for disease genes.
  numerator <- sum(influence_values[valid_disease_genes])
  
  # Calculate denominator: sum of influences for non-disease genes.
  denominator <- sum(influence_values[non_disease_genes])
  
  if (denominator == 0) {
    return(0)
  }
  
  return(numerator / denominator)
}

# Header: first column is Tissue-Name, then one column per candidate drug.
efficacy_matrix <- list()
efficacy_matrix[[1]] <- c("Tissue-Name", candidate_drugs)

row_index <- 2  # Row 1 is the header.
for (mapped_file in mapped_files) {
  filename    <- basename(mapped_file)  # e.g., "mapped_brain_cortex_20.csv"
  tissue_name <- sub(paste0("^mapped_(.*)_", age_group, "\\.csv$"), "\\1", filename)
  
  cat(sprintf("\nProcessing tissue: %s (file=%s)\n", tissue_name, filename))
  
  edges <- read.csv(mapped_file, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(edges) == 0) {
    # No edges: assign efficacy = 0 for all candidate drugs.
    row_values <- c(tissue_name, rep("0", length(candidate_drugs)))
    efficacy_matrix[[row_index]] <- row_values
    row_index <- row_index + 1
    cat("  --> No edges found, all efficacy=0.\n")
    next
  }
  
  # Build the tissue network graph.
  g <- graph_from_data_frame(edges, directed = FALSE)
  all_genes <- V(g)$name
  
  # non_disease_genes: all genes in the network not in the full disease gene list.
  non_disease_genes <- setdiff(all_genes, disease_genes)
  
  # For each candidate drug, compute the efficacy.
  drug_vals <- character(length(candidate_drugs))
  for (i in seq_along(candidate_drugs)) {
    drug <- candidate_drugs[i]
    targets <- candidate_drug_targets_list[[drug]]
    
    # Check if at least one target gene is present in this tissue network.
    if (length(intersect(targets, all_genes)) == 0) {
      cat(sprintf("  --> Warning: none of the targets for '%s' found in tissue %s. Setting efficacy=0.\n",
                  drug, tissue_name))
      drug_vals[i] <- "0"
    } else {
      eff_val <- calc_efficacy(g, targets, disease_genes, non_disease_genes, lambda)
      drug_vals[i] <- as.character(eff_val)
    }
  }
  
  # Build the row for this tissue and add it to the results.
  row_values <- c(tissue_name, drug_vals)
  efficacy_matrix[[row_index]] <- row_values
  row_index <- row_index + 1
}

cat(sprintf("\nWriting CSV to: %s\n", output_file))
file_conn <- file(output_file, open = "w")
for (row_vec in efficacy_matrix) {
  line_str <- paste(row_vec, collapse = ",")
  writeLines(line_str, con = file_conn)
}
close(file_conn)
cat("Done.\n")
