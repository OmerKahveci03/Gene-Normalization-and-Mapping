# drug_efficacy.R - takes roughly 9 minutes
# -----------------------------------------------------------
# Input:
#   - Disease genes S (read from filtered_neurodegenerative_diseases.csv) for Alzheimer Disease.
#   - Drug candidate information from filtered_drug_gene_dataset.xlsx:
#         * "Sorted" sheet: Column A = drug name, Column B = gene count (already sorted descending).
#         * "Raw" sheet: Column A = drug name, Column B = gene name (there may be duplicates).
#
#   For each candidate drug (from Sorted with gene count >= 20), we retrieve its target genes (from Raw),
#   compute the intersection with the Alzheimer Disease genes, and select the drug with the highest intersection.
#   Then we pick a representative drug target gene from that drug (preferring one that is in the disease set)
#   and use that for efficacy calculation.
#
#   The efficacy formula for a drug target gene d:
#      efficacy(d, S) = (sum_{s in S} influence_g1_g2(d, s)) / (sum_{n in N} influence_g1_g2(d, n))
#      where N = (all genes in the tissue network) \ S.
#
#   Output:
#     - A CSV file (data/drug/efficacy.csv) with tissues in rows and one column for the selected drug.
#
library(igraph)
library(readxl)
source("utility.R")

# ------------------------------------------------------------------------------
# 1. Read Disease Genes for 'Alzheimer Disease'
# ------------------------------------------------------------------------------
nd_file <- file.path("..", "data", "drug", "filtered_neurodegenerative_diseases.csv")
if (!file.exists(nd_file)) {
  stop(sprintf("Neurodegenerative diseases file not found: %s", nd_file))
}

# Assuming the file has no header; if it does, set header = TRUE and adjust column names.
nd_data <- read.csv(nd_file, header = FALSE, stringsAsFactors = FALSE)
# Column 1: gene name, Column 2: disease name.
disease_genes <- nd_data[nd_data[, 2] == "Alzheimer Disease", 1]
if (length(disease_genes) == 0) {
  stop("No genes found for 'Alzheimer Disease' in the neurodegenerative diseases file.")
}

# (Optional) Clean the gene names (trim whitespace and/or convert case if needed)
disease_genes <- trimws(disease_genes)

# ------------------------------------------------------------------------------
# 2. Read Drug Candidate Information from Excel
# ------------------------------------------------------------------------------
drug_excel_file <- file.path("..", "data", "drug", "filtered_drug_gene_dataset.xlsx")
if (!file.exists(drug_excel_file)) {
  stop(sprintf("Drug Excel file not found: %s", drug_excel_file))
}

# Read the "Sorted" sheet (Column 1: drug name, Column 2: gene count).
sorted_data <- read_excel(drug_excel_file, sheet = "Sorted", col_names = FALSE)
# Only consider drugs with gene count >= 20.
candidate_sorted <- sorted_data[sorted_data[[2]] >= 20, ]
if(nrow(candidate_sorted) == 0) {
  stop("No candidate drugs with gene count >= 20 found in the Sorted sheet.")
}

# Read the "Raw" sheet (Column 1: drug name, Column 2: gene name).
raw_data <- read_excel(drug_excel_file, sheet = "Raw", col_names = FALSE)

# For each candidate drug, compute the intersection size between its target genes and disease_genes.
candidate_intersections <- sapply(candidate_sorted[[1]], function(drug) {
  # Retrieve target genes for this drug from the Raw sheet.
  targets <- raw_data[[2]][raw_data[[1]] == drug]
  length(intersect(targets, disease_genes))
})

# Select the candidate drug with the highest intersection size.
max_intersection <- max(candidate_intersections)
selected_index <- which(candidate_intersections == max_intersection)[1]
selected_drug <- candidate_sorted[[1]][selected_index]
# Retrieve all target genes for the selected drug.
selected_drug_targets <- raw_data[[2]][raw_data[[1]] == selected_drug]
# Compute the intersection with disease_genes.
selected_drug_intersection <- intersect(selected_drug_targets, disease_genes)
# Pick a representative drug target gene.
if (length(selected_drug_intersection) > 0) {
  selected_drug_gene <- selected_drug_intersection[1]
} else {
  selected_drug_gene <- selected_drug_targets[1]
}

cat(sprintf("Selected drug: %s (using target gene: %s) with intersection count: %d\n",
            selected_drug, selected_drug_gene, max_intersection))

# ------------------------------------------------------------------------------
# 3. Configurable Variables for the Tissue Networks
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 4. Prepare the Results Container
# ------------------------------------------------------------------------------
# The header now has "Tissue-Name" and the selected drug name.
efficacy_matrix <- list()
efficacy_matrix[[1]] <- c("Tissue-Name", selected_drug)

# ------------------------------------------------------------------------------
# 5. Updated Helper Function: Calculate Efficacy Using a Set of Drug Target Genes
# ------------------------------------------------------------------------------
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
  
  numerator <- sum(sapply(valid_disease_genes, function(s) {
    influence_S_g(graph, valid_drug_targets, s, lambda)
  }))
  
  denominator <- sum(sapply(non_disease_genes, function(n) {
    influence_S_g(graph, valid_drug_targets, n, lambda)
  }))
  
  if (denominator == 0) {
    return(0)
  }
  
  return(numerator / denominator)
}

# ------------------------------------------------------------------------------
# 6. Loop Over Each Tissue CSV and Compute the Efficacy for the Selected Drug Targets
# ------------------------------------------------------------------------------
row_index <- 2  # first row is the header
for (mapped_file in mapped_files) {
  filename    <- basename(mapped_file)  # e.g., "mapped_brain_cortex_20.csv"
  tissue_name <- sub(paste0("^mapped_(.*)_", age_group, "\\.csv$"), "\\1", filename)
  
  cat(sprintf("\nProcessing tissue: %s (file=%s)\n", tissue_name, filename))
  
  edges <- read.csv(mapped_file, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(edges) == 0) {
    row_values <- c(tissue_name, "0")
    efficacy_matrix[[row_index]] <- row_values
    row_index <- row_index + 1
    cat("  --> No edges found, efficacy=0.\n")
    next
  }
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  all_genes <- V(g)$name
  
  # non_disease_genes: all genes in the network not in the full disease list.
  non_disease_genes <- setdiff(all_genes, disease_genes)
  
  if (length(intersect(selected_drug_targets, all_genes)) == 0) {
    cat(sprintf("  --> Warning: none of the targets for '%s' found in tissue %s. Setting efficacy=0.\n", 
                selected_drug, tissue_name))
    efficacy_val <- 0
  } else {
    efficacy_val <- calc_efficacy(g, selected_drug_targets, disease_genes, non_disease_genes, lambda)
  }
  
  row_values <- c(tissue_name, as.character(efficacy_val))
  efficacy_matrix[[row_index]] <- row_values
  row_index <- row_index + 1
}


# ------------------------------------------------------------------------------
# 7. Write the Efficacy Matrix to CSV
# ------------------------------------------------------------------------------
cat(sprintf("\nWriting CSV to: %s\n", output_file))
file_conn <- file(output_file, open = "w")
for (row_vec in efficacy_matrix) {
  line_str <- paste(row_vec, collapse = ",")
  writeLines(line_str, con = file_conn)
}
close(file_conn)
cat("Done.\n")
