# drug_efficacy.R
# -----------------------------------------------------------
# Input:
#   - set of disease genes S = {s1, s2, s3,... sj}
#   - set of drug target genes D = {d1, d2, d3,... di}
#   - age_group
#   - for each tissue within data/mapped/{age_group}, read the network CSV
#
#   The "influence_g1_g2" function is sourced from utility.R
#   (Make sure utility.R has influence_g1_g2 defined)
#
# Efficacy formula for a single drug di:
#   numerator   = sum(inf(di, s_j))  over all s_j in S
#   denominator = sum(inf(di, n_j))  over all n_j in N = (V \ S)
#   efficacy(di, S) = numerator / denominator
#
# Output:
#   - A file data/drug/efficacy.csv listing tissues in rows and drug efficacies in columns.
#     The first row is: "Tissue Name,AATF,ABL1,AES,AHR" (for example).

library(igraph)
source("utility.R")

# Configurable Variables
disease_genes     <- c("CDH2", "CXCR4")              # Example disease genes
drug_target_genes <- c("AATF", "ABL1", "AES", "AHR") # Example drug target genes
age_group         <- 20
lambda            <- 1

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

# Stores the results
efficacy_matrix <- list()
efficacy_matrix[[1]] <- c("Tissue-Name", drug_target_genes)

# Helper function: calculates efficacy for a single drug w.r.t. disease genes
calc_efficacy <- function(graph, d_gene, disease_genes, non_disease_genes, lambda = 1) {
  # numerator = sum of inf(d_gene, disease_gene)
  numerator <- sum(sapply(disease_genes, function(s) {
    influence_g1_g2(graph, d_gene, s, lambda)
  }))
  
  # denominator = sum of inf(d_gene, non_disease_genes)
  denominator <- sum(sapply(non_disease_genes, function(n) {
    influence_g1_g2(graph, d_gene, n, lambda)
  }))
  
  # If denominator=0 or missing => efficacy=0
  if (denominator == 0) {
    return(0)
  }
  return(numerator / denominator)
}

# ------------------------------------------------------------------------------
# Loop over each CSV in "mapped_dir" and compute a new row for each tissue
# ------------------------------------------------------------------------------

row_index <- 2  # Because row 1 is the header
for (mapped_file in mapped_files) {
  filename     <- basename(mapped_file)  # e.g. "mapped_brain_cortex_20.csv"
  tissue_name  <- sub(paste0("^mapped_(.*)_", age_group, "\\.csv$"), "\\1", filename)
  
  cat(sprintf("\nProcessing tissue: %s (file=%s)\n", tissue_name, filename))
  
  edges <- read.csv(mapped_file, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(edges) == 0) {
    # No edges => all drug efficacy = 0
    row_values <- c(tissue_name, rep("0", length(drug_target_genes)))
    efficacy_matrix[[row_index]] <- row_values
    row_index <- row_index + 1
    cat(sprintf("  --> No edges found, all efficacy=0.\n"))
    next
  }
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  all_genes <- V(g)$name
  
  # Non-disease genes = (all_genes - disease_genes)
  non_disease_genes <- setdiff(all_genes, disease_genes)
  
  # For each drug, compute efficacy
  drug_vals <- character(length(drug_target_genes))
  for (i in seq_along(drug_target_genes)) {
    d_gene <- drug_target_genes[i]
    if (!d_gene %in% all_genes) {
      # If drug gene doesn't exist in this network => 0
      cat(sprintf("  --> Warning: drug gene '%s' not found. Setting efficacy=0\n", d_gene))
      drug_vals[i] <- "0"
    } else {
      eff_val <- calc_efficacy(g, d_gene, disease_genes, non_disease_genes, lambda)
      drug_vals[i] <- as.character(eff_val)
    }
  }
  
  # Build the row for this tissue
  row_values <- c(tissue_name, drug_vals)
  efficacy_matrix[[row_index]] <- row_values
  row_index <- row_index + 1
}

# ------------------------------------------------------------------------------
# Write the "efficacy_matrix" to CSV manually
# ------------------------------------------------------------------------------

cat(sprintf("\nWriting CSV to: %s\n", output_file))
file_conn <- file(output_file, open = "w")

for (row_vec in efficacy_matrix) {
  line_str <- paste(row_vec, collapse = ",")
  writeLines(line_str, con = file_conn)
}

close(file_conn)
cat("Done.\n")
