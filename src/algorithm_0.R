# algorithm_0.R

library(igraph)
source("utility.R")  # Make sure compute_influence_list() & age_based_influence() are in utility.R

# -----------------------------
# Configurable Values
# -----------------------------
tissue_name  <- "brain_cortex"
lambda       <- 1
epsilon      <- 0.005

age_groups <- c(20, 30, 40, 50, 60, 70)

# Starting gene set
gene_set <- character(0)

all_genes <- character(0)

# get all_gets
for (ag in age_groups) {
  mapping_file <- file.path(
    "..", "data", "mapped",
    as.character(ag),
    paste0("mapped_", tissue_name, "_", ag, ".csv")
  )
  
  if (!file.exists(mapping_file)) {
    cat(sprintf("Age group %d: File does not exist => %s\n", ag, mapping_file))
    next
  }
  
  edges <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(edges) == 0) {
    cat(sprintf("Age group %d: No edges found.\n", ag))
    next
  }
  genes_in_file <- unique(c(edges[[1]], edges[[2]]))
  all_genes     <- c(all_genes, genes_in_file)
}

all_genes <- unique(all_genes)
cat(sprintf("Total unique genes across mapped files: %d\n", length(all_genes)))

# Filter out any genes already in gene_set (though it's empty now)
candidate_genes <- setdiff(all_genes, gene_set)

# ------------------------------------------------------------
# 2) One iteration of "greedy": test each candidate gene
# ------------------------------------------------------------
score_dict <- list()  # store key-value pairs: gene -> influence_score

cat("\nStarting single iteration of greedy selection...\n")
start_time <- Sys.time()

for (idx in seq_along(candidate_genes)) {
  g <- candidate_genes[idx]
  
  # Temporarily add gene g to gene_set
  temp_set <- c(gene_set, g)
  
  # Compute the influence_list for this temp_set
  influence_list <- compute_influence_list(
    age_groups  = age_groups,
    tissue_name = tissue_name,
    gene_set    = temp_set,
    lambda      = lambda,
    data_dir    = ".."
  )
  
  # Get the final age-based influence metric
  ab_infl_val <- age_based_influence(influence_list, epsilon)
  
  # Store in our dictionary
  score_dict[[g]] <- ab_infl_val
  
  # Print progress
  cat(sprintf("Gene [%d/%d] completed.\n", idx, length(candidate_genes)))
}

end_time <- Sys.time()
time_diff_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

# ------------------------------------------------------------
# 3) Pick the gene with the highest age-based influence value
# ------------------------------------------------------------
best_gene <- NA
best_val  <- -Inf

for (gene_name in names(score_dict)) {
  val <- score_dict[[gene_name]]
  if (!is.na(val) && val > best_val) {
    best_val  <- val
    best_gene <- gene_name
  }
}

# ------------------------------------------------------------
# 4) Permanently add that best_gene to gene_set
# ------------------------------------------------------------
if (!is.na(best_gene)) {
  gene_set <- c(gene_set, best_gene)
}

# ------------------------------------------------------------
# 5) Clear the dictionary
# ------------------------------------------------------------
score_dict <- list()

# ------------------------------------------------------------
# 6) Print the chosen gene, its influence val, and timing
# ------------------------------------------------------------
cat(sprintf("\nGreedy pick:\n"))
cat(sprintf("  Chosen Gene: %s\n", best_gene))
cat(sprintf("  Influence:   %f\n", best_val))
cat(sprintf("  Time (secs): %.3f\n", time_diff_sec))
cat(sprintf("  New gene_set size = %d\n", length(gene_set)))
cat(sprintf("  Genes now in gene_set: %s\n", paste(gene_set, collapse=", ")))
