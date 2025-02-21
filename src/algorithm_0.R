library(igraph)

# -----------------------------
# Configurable Values
# -----------------------------
tissue_name  <- "brain_cortex"
lambda       <- 1
epsilon      <- 0.005
age_groups   <- c(20, 30, 40, 50, 60, 70)
# Instead of a fixed number of iterations, we loop until no candidate improves the metric.

# -----------------------------
# Start overall timer
# -----------------------------
start_time_total <- Sys.time()

# -----------------------------
# Load TRRUST file to get all unique gene symbols
# -----------------------------
trrust_file <- file.path("..", "data", "raw", "trrust_rawdata.human.tsv")
if (!file.exists(trrust_file)) {
  stop("TRRUST file not found at: ", trrust_file)
}
# Read the TRRUST file (assumes tab-separated; first two columns are gene symbols)
trrust <- read.delim(trrust_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Get all unique genes from both columns
trrust_genes <- unique(c(trrust[[1]], trrust[[2]]))

# -----------------------------
# Pre-load Graphs and Collect Source Genes from Mapped Files
# -----------------------------
graphs <- list()
mapped_source_genes <- character(0)  # union of all source genes (column A) from mapped files

for (ag in age_groups) {
  mapping_file <- file.path("..", "data", "mapped", as.character(ag),
                            paste0("mapped_", tissue_name, "_", ag, ".csv"))
  
  if (!file.exists(mapping_file)) {
    cat(sprintf("Age group %d: File does not exist => %s\n", ag, mapping_file))
    graphs[[as.character(ag)]] <- NULL
    next
  }
  
  edges <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(edges) == 0) {
    cat(sprintf("Age group %d: No edges found.\n", ag))
    graphs[[as.character(ag)]] <- NULL
    next
  }
  
  # Build the directed graph (column A influences column B)
  g <- graph_from_data_frame(edges, directed = TRUE)
  graphs[[as.character(ag)]] <- g
  
  # Update candidate source genes: only genes from column A (influencers)
  source_genes <- unique(edges[[1]])
  mapped_source_genes <- union(mapped_source_genes, source_genes)
}

cat(sprintf("Total unique source genes from mapped files: %d\n", length(mapped_source_genes)))

# Create final candidate gene set: only genes from TRRUST that are found as influencers (column A)
all_candidate_genes <- intersect(trrust_genes, mapped_source_genes)
cat(sprintf("Total candidate genes (from TRRUST and mapped sources): %d\n", length(all_candidate_genes)))

# -----------------------------
# Influence Computation Functions
# -----------------------------
# Computes the average influence of set S on all nodes in graph g.
# Influence(S, G) = (sum_{g in G} e^(-lambda * min_{s in S} dist(s, g))) / |G|
compute_influence_preloaded <- function(g, S, lambda = 1) {
  if (length(S) == 0) return(0)
  
  all_nodes <- V(g)$name
  # Filter S to only include valid vertices in the graph.
  valid_S <- intersect(S, all_nodes)
  if (length(valid_S) == 0) return(0)
  
  # Compute the full distance matrix from valid_S to all nodes in g.
  dist_mat <- distances(g, v = valid_S, to = all_nodes)
  # For each node g in G, take the minimum distance from any s in valid_S.
  min_dists <- apply(dist_mat, 2, min)
  influence_values <- exp(-lambda * min_dists)
  
  return(mean(influence_values))
}

# Computes the influence list (one per age group) across all pre-loaded graphs.
compute_influence_list_preloaded <- function(graphs, gene_set, lambda = 1) {
  sapply(graphs, function(g) {
    if (is.null(g)) return(NA)
    compute_influence_preloaded(g, gene_set, lambda)
  })
}

# Computes the age-based influence from the influence list.
age_based_influence <- function(influence_list, epsilon = 0.005) {
  val <- 0
  # Compare age groups 2 through 5 with later groups.
  for (i in 2:5) {
    for (j in (i+1):6) {
      if (is.na(influence_list[i]) || is.na(influence_list[j])) {
        next
      }
      if (influence_list[i] + epsilon < influence_list[j]) {
        val <- val + 1
      } else if (influence_list[j] + epsilon < influence_list[i]) {
        val <- val - 1
      }
    }
  }
  return(val / 15)
}

# -----------------------------
# Greedy Search for Maximum Age-Based Influence
# -----------------------------
cat("\n=== Greedy Search for Maximum Age-Based Influence ===\n")
gene_set_max <- character(0)
current_metric_max <- age_based_influence(compute_influence_list_preloaded(graphs, gene_set_max, lambda), epsilon)
improved <- TRUE
iter_max <- 0

while (improved) {
  iter_max <- iter_max + 1
  cat(sprintf("\nMax Iteration %d\n", iter_max))
  candidate_genes <- setdiff(all_candidate_genes, gene_set_max)
  best_candidate <- NA
  best_candidate_metric <- current_metric_max
  
  for (candidate in candidate_genes) {
    temp_set <- c(gene_set_max, candidate)
    influence_list <- compute_influence_list_preloaded(graphs, temp_set, lambda)
    new_metric <- age_based_influence(influence_list, epsilon)
    
    if (new_metric > best_candidate_metric) {
      best_candidate_metric <- new_metric
      best_candidate <- candidate
    }
  }
  
  if (!is.na(best_candidate) && best_candidate_metric > current_metric_max) {
    gene_set_max <- c(gene_set_max, best_candidate)
    current_metric_max <- best_candidate_metric
    cat(sprintf("  Chosen Gene: %s  New Metric: %f\n", best_candidate, current_metric_max))
    cat(sprintf("  Gene Set: %s\n", paste(gene_set_max, collapse = ", ")))
  } else {
    cat("No candidate gene improves the metric further (maximization).\n")
    improved <- FALSE
  }
}

# -----------------------------
# Greedy Search for Minimum Age-Based Influence (Most Negative)
# -----------------------------
cat("\n=== Greedy Search for Minimum Age-Based Influence ===\n")
gene_set_min <- character(0)
current_metric_min <- age_based_influence(compute_influence_list_preloaded(graphs, gene_set_min, lambda), epsilon)
improved <- TRUE
iter_min <- 0

while (improved) {
  iter_min <- iter_min + 1
  cat(sprintf("\nMin Iteration %d\n", iter_min))
  candidate_genes <- setdiff(all_candidate_genes, gene_set_min)
  best_candidate <- NA
  best_candidate_metric <- current_metric_min
  
  for (candidate in candidate_genes) {
    temp_set <- c(gene_set_min, candidate)
    influence_list <- compute_influence_list_preloaded(graphs, temp_set, lambda)
    new_metric <- age_based_influence(influence_list, epsilon)
    
    if (new_metric < best_candidate_metric) {  # looking for more negative value
      best_candidate_metric <- new_metric
      best_candidate <- candidate
    }
  }
  
  if (!is.na(best_candidate) && best_candidate_metric < current_metric_min) {
    gene_set_min <- c(gene_set_min, best_candidate)
    current_metric_min <- best_candidate_metric
    cat(sprintf("  Chosen Gene: %s  New Metric: %f\n", best_candidate, current_metric_min))
    cat(sprintf("  Gene Set: %s\n", paste(gene_set_min, collapse = ", ")))
  } else {
    cat("No candidate gene improves the metric further (minimization).\n")
    improved <- FALSE
  }
}

# -----------------------------
# Plotting the Final Influence Values for Both Gene Sets
# -----------------------------
final_influence_list_max <- compute_influence_list_preloaded(graphs, gene_set_max, lambda)
final_influence_list_min <- compute_influence_list_preloaded(graphs, gene_set_min, lambda)

# Plot for Maximization
plot(age_groups, final_influence_list_max, type = "b", pch = 19, col = "blue",
     xlab = "Age Group", ylab = "Influence(S, V)",
     main = "Influence by Age Group (Maximized Influence)",
     ylim = c(0, max(final_influence_list_max, na.rm = TRUE) * 1.1))
# Plot for Minimization
plot(age_groups, final_influence_list_min, type = "b", pch = 19, col = "red",
     xlab = "Age Group", ylab = "Influence(S, V)",
     main = "Influence by Age Group (Minimized Influence)",
     ylim = c(0, max(final_influence_list_min, na.rm = TRUE) * 1.1))

# -----------------------------
# Report Final Results and Total Runtime
# -----------------------------
cat("\n=== Final Results ===\n")
cat(sprintf("Maximum Influence Gene Set: %s\n", paste(gene_set_max, collapse = ", ")))
cat(sprintf("Final Maximum Age-Based Influence: %f\n", current_metric_max))
cat(sprintf("Minimum Influence Gene Set: %s\n", paste(gene_set_min, collapse = ", ")))
cat(sprintf("Final Minimum Age-Based Influence: %f\n", current_metric_min))

end_time_total <- Sys.time()
total_time_sec <- as.numeric(difftime(end_time_total, start_time_total, units = "secs"))
cat(sprintf("\nTotal runtime: %.3f seconds\n", total_time_sec))
