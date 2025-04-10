# Total runtime: 1738.537 seconds
library(igraph)
library(writexl)

# -----------------------------
# Determine Tissue List from mapped files for age 20
# -----------------------------
tissue_dir <- file.path("..", "data", "mapped", "20")
mapped_files <- list.files(path = tissue_dir, pattern = "^mapped_.*_20\\.csv$", full.names = FALSE)
# Extract tissue names from file names of the form "mapped_{tissue}_20.csv"
tissues <- unique(gsub("^mapped_(.*)_20\\.csv$", "\\1", mapped_files))
cat(sprintf("Found %d tissues: %s\n", length(tissues), paste(tissues, collapse = ", ")))

# -----------------------------
# Configurable Values
# -----------------------------
lambda <- 1
epsilon <- 0.005
age_groups <- c(20, 30, 40, 50, 60, 70)
s_lambda <- 5  # Constant to raise improvement scores

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
trrust <- read.delim(trrust_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
trrust_genes <- unique(c(trrust[[1]], trrust[[2]]))

# -----------------------------
# Define Incremental Update and stochastic Search Functions with Stochasticity
# -----------------------------
initialize_distances <- function(g) {
  if (is.null(g)) return(NULL)
  rep(Inf, length(V(g)$name))
}

# Updated: If candidate is not in the graph, leave distances unchanged.
update_influence_for_graph <- function(g, d_current, candidate, lambda) {
  all_nodes <- V(g)$name
  if (!(candidate %in% all_nodes)) {
    # Candidate not in the graph: return current distances and influence unchanged.
    return(list(d_new = d_current, influence_new = mean(exp(-lambda * d_current))))
  }
  d_candidate <- as.vector(distances(g, v = candidate, to = all_nodes))
  d_new <- pmin(d_current, d_candidate)
  influence_new <- mean(exp(-lambda * d_new))
  list(d_new = d_new, influence_new = influence_new)
}

age_based_influence <- function(influence_list, epsilon = 0.005) {
  val <- 0
  for (i in 2:5) {
    for (j in (i+1):6) {
      if (is.na(influence_list[i]) || is.na(influence_list[j])) next
      if (influence_list[i] + epsilon < influence_list[j]) {
        val <- val + 1
      } else if (influence_list[j] + epsilon < influence_list[i]) {
        val <- val - 1
      }
    }
  }
  return(val / 15)
}

stochastic_search <- function(mode = c("max", "min"), candidate_genes, graphs, lambda, epsilon, s_lambda) {
  mode <- match.arg(mode)
  gene_set <- character(0)
  current_d <- lapply(graphs, initialize_distances)
  current_influences <- mapply(function(g, d) {
    if (is.null(g)) return(NA)
    mean(exp(-lambda * d))
  }, graphs, current_d)
  current_metric <- age_based_influence(current_influences, epsilon)
  
  improved <- TRUE
  iteration <- 0
  
  while (improved) {
    iteration <- iteration + 1
    cat(sprintf("\n%s Iteration %d\n", toupper(mode), iteration))
    remaining_candidates <- setdiff(candidate_genes, gene_set)
    
    # Storage for candidates that improve the metric.
    candidate_list <- list()
    candidate_updates_list <- list()
    candidate_weights <- numeric()
    
    for (candidate in remaining_candidates) {
      new_influences <- numeric(length(graphs))
      candidate_updates <- vector("list", length(graphs))
      
      for (i in seq_along(graphs)) {
        g <- graphs[[i]]
        if (is.null(g)) {
          new_influences[i] <- NA
          candidate_updates[[i]] <- NULL
          next
        }
        res <- update_influence_for_graph(g, current_d[[i]], candidate, lambda)
        candidate_updates[[i]] <- res$d_new
        new_influences[i] <- res$influence_new
      }
      new_metric <- age_based_influence(new_influences, epsilon)
      
      # Only consider candidates that improve the metric.
      if ((mode == "max" && new_metric > current_metric) ||
          (mode == "min" && new_metric < current_metric)) {
        
        # When current_metric is 0, use absolute difference (only meaningful for max mode)
        if (current_metric == 0) {
          if (mode == "max") {
            improvement_factor <- new_metric  # absolute improvement
            weight <- improvement_factor * s_lambda  # multiply by s_lambda
          } else {
            # For min mode, 0 is the optimal and cannot be improved upon.
            next
          }
        } else {
          # Otherwise, use percentage change.
          if (mode == "max") {
            improvement_factor <- new_metric / current_metric
          } else {
            improvement_factor <- current_metric / new_metric
          }
          weight <- improvement_factor^s_lambda
        }
        
        # Only store candidate if the computed values are valid.
        if (is.na(improvement_factor) || !is.finite(improvement_factor) ||
            weight <= 0 || is.na(weight) || !is.finite(weight)) {
          next
        }
        
        candidate_list <- c(candidate_list, candidate)
        candidate_updates_list[[length(candidate_list)]] <- candidate_updates
        candidate_weights <- c(candidate_weights, weight)
      }
    }
    
    if (length(candidate_list) == 0 || sum(candidate_weights > 0) == 0) {
      cat(sprintf("No candidate gene improves the metric further (%s).\n", mode))
      improved <- FALSE
      break
    }
    
    # Randomly select one candidate using the computed weights.
    selected_index <- sample(seq_along(candidate_list), size = 1, prob = candidate_weights)
    best_candidate <- candidate_list[[selected_index]]
    best_candidate_updates <- candidate_updates_list[[selected_index]]
    
    # Update the current influences using the selected candidate.
    new_influences <- numeric(length(graphs))
    for (i in seq_along(graphs)) {
      g <- graphs[[i]]
      if (is.null(g)) {
        new_influences[i] <- NA
        next
      }
      res <- update_influence_for_graph(g, current_d[[i]], best_candidate, lambda)
      new_influences[i] <- res$influence_new
    }
    new_metric <- age_based_influence(new_influences, epsilon)
    
    # Update the gene set and current distances.
    gene_set <- c(gene_set, best_candidate)
    current_metric <- new_metric
    for (i in seq_along(graphs)) {
      if (!is.null(graphs[[i]]))
        current_d[[i]] <- best_candidate_updates[[i]]
    }
    cat(sprintf("  Chosen Gene: %s  New Metric: %f\n", best_candidate, current_metric))
    cat(sprintf("  Gene Set: %s\n", paste(gene_set, collapse = ", ")))
  }
  
  final_influences <- mapply(function(g, d) {
    if (is.null(g)) return(NA)
    mean(exp(-lambda * d))
  }, graphs, current_d)
  list(gene_set = gene_set, metric = current_metric, influences = final_influences)
}

# -----------------------------
# Prepare to store Excel results for each tissue
# -----------------------------
max_results_list <- list()
min_results_list <- list()

# -----------------------------
# Loop over all tissues
# -----------------------------
for (tissue in tissues) {
  cat(sprintf("\n=============================\nProcessing Tissue: %s\n", tissue))
  
  # For this tissue, load graphs for each age group
  graphs_tissue <- list()
  mapped_source_genes_tissue <- character(0)
  
  for (ag in age_groups) {
    mapping_file <- file.path("..", "data", "mapped", as.character(ag),
                              paste0("mapped_", tissue, "_", ag, ".csv"))
    if (!file.exists(mapping_file)) {
      cat(sprintf("Tissue %s, Age %d: File does not exist => %s\n", tissue, ag, mapping_file))
      graphs_tissue[[as.character(ag)]] <- NULL
      next
    }
    edges <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)
    if (nrow(edges) == 0) {
      cat(sprintf("Tissue %s, Age %d: No edges found.\n", tissue, ag))
      graphs_tissue[[as.character(ag)]] <- NULL
      next
    }
    g <- graph_from_data_frame(edges, directed = TRUE)
    graphs_tissue[[as.character(ag)]] <- g
    source_genes <- unique(edges[[1]])
    mapped_source_genes_tissue <- union(mapped_source_genes_tissue, source_genes)
  }
  
  cat(sprintf("Tissue %s: Total unique source genes from mapped files: %d\n",
              tissue, length(mapped_source_genes_tissue)))
  
  candidate_genes_tissue <- intersect(trrust_genes, mapped_source_genes_tissue)
  cat(sprintf("Tissue %s: Total candidate genes: %d\n", tissue, length(candidate_genes_tissue)))
  
  if (length(candidate_genes_tissue) == 0) {
    cat(sprintf("Tissue %s: No candidate genes available. Skipping.\n", tissue))
    next
  }
  
  cat(sprintf("Tissue %s: Running stochastic search for MAX influence...\n", tissue))
  result_max <- stochastic_search("max", candidate_genes_tissue, graphs_tissue, lambda, epsilon, s_lambda)
  cat(sprintf("Tissue %s: Running stochastic search for MIN influence...\n", tissue))
  result_min <- stochastic_search("min", candidate_genes_tissue, graphs_tissue, lambda, epsilon, s_lambda)
  
  output_folder <- file.path("..", "data", "age_indexes")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  max_plot_file <- file.path(output_folder, sprintf("max_%s_plot.png", tissue))
  png(filename = max_plot_file, width = 800, height = 600)
  plot(age_groups, result_max$influences, type = "b", pch = 19, col = "blue",
       xlab = "Age Group", ylab = "Influence(S, V)",
       main = sprintf("Influence by Age Group (Maximized) - %s", tissue),
       ylim = c(0, max(result_max$influences, na.rm = TRUE) * 1.1))
  dev.off()
  
  min_plot_file <- file.path(output_folder, sprintf("min_%s_plot.png", tissue))
  png(filename = min_plot_file, width = 800, height = 600)
  plot(age_groups, result_min$influences, type = "b", pch = 19, col = "red",
       xlab = "Age Group", ylab = "Influence(S, V)",
       main = sprintf("Influence by Age Group (Minimized) - %s", tissue),
       ylim = c(0, max(result_min$influences, na.rm = TRUE) * 1.1))
  dev.off()
  
  max_row <- c(tissue, result_max$metric, result_max$gene_set)
  min_row <- c(tissue, result_min$metric, result_min$gene_set)
  
  max_results_list[[length(max_results_list) + 1]] <- max_row
  min_results_list[[length(min_results_list) + 1]] <- min_row
  
  cat(sprintf("Tissue %s processed. MAX metric: %f, MIN metric: %f\n", tissue,
              result_max$metric, result_min$metric))
}

fill_row <- function(x, target_length) {
  length(x) <- target_length
  x
}

max_ncols <- max(sapply(max_results_list, length))
min_ncols <- max(sapply(min_results_list, length))

max_mat <- do.call(rbind, lapply(max_results_list, fill_row, target_length = max_ncols))
min_mat <- do.call(rbind, lapply(min_results_list, fill_row, target_length = min_ncols))

max_df <- as.data.frame(max_mat, stringsAsFactors = FALSE)
min_df <- as.data.frame(min_mat, stringsAsFactors = FALSE)

max_excel_file <- file.path("..", "data", "age_indexes", "a1_max_age_indexes.xlsx")
min_excel_file <- file.path("..", "data", "age_indexes", "a1_min_age_indexes.xlsx")

write_xlsx(max_df, path = max_excel_file, col_names = FALSE)
write_xlsx(min_df, path = min_excel_file, col_names = FALSE)

end_time_total <- Sys.time()
total_time_sec <- as.numeric(difftime(end_time_total, start_time_total, units = "secs"))
cat(sprintf("\nTotal runtime: %.3f seconds\n", total_time_sec))
