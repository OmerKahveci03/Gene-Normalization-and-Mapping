# Total runtime: 9520.918 seconds
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
s_lambda <- 5  # Constant to scale improvement scores
num_runs <- 5  # Number of stochastic runs per tissue

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
# Define Incremental Update and Stochastic Search Functions
# -----------------------------
initialize_distances <- function(g) {
  if (is.null(g)) return(NULL)
  rep(Inf, length(V(g)$name))
}

update_influence_for_graph <- function(g, d_current, candidate, lambda) {
  all_nodes <- V(g)$name
  if (!(candidate %in% all_nodes)) {
    # Candidate not in the graph: distances and influence unchanged.
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
      
      if ((mode == "max" && new_metric > current_metric) ||
          (mode == "min" && new_metric < current_metric)) {
        
        if (current_metric == 0) {
          if (mode == "max") {
            improvement_factor <- new_metric  # absolute improvement
            weight <- improvement_factor * s_lambda  # multiply by s_lambda
          } else {
            # In min mode, 0 is optimal; skip.
            next
          }
        } else {
          if (mode == "max") {
            improvement_factor <- new_metric / current_metric
          } else {
            improvement_factor <- current_metric / new_metric
          }
          weight <- improvement_factor^s_lambda
        }
        
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
    
    selected_index <- sample(seq_along(candidate_list), size = 1, prob = candidate_weights)
    best_candidate <- candidate_list[[selected_index]]
    best_candidate_updates <- candidate_updates_list[[selected_index]]
    
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
# Prepare lists to store Excel results for each mode and run
# Each will be a list of length = num_runs where each element is a list of rows (one per tissue)
# Also, for each tissue we will store the best row (based on metric) among runs.
# -----------------------------
max_runs_list <- vector("list", num_runs)
min_runs_list <- vector("list", num_runs)
for (i in 1:num_runs) {
  max_runs_list[[i]] <- list()
  min_runs_list[[i]] <- list()
}
max_best_list <- list()
min_best_list <- list()

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
  
  ## For each tissue, run stochastic search num_runs times for MAX and MIN modes.
  best_max_metric <- -Inf   # For max mode, higher is better.
  best_min_metric <- Inf    # For min mode, lower is better.
  best_max_row <- NULL
  best_min_row <- NULL
  
  for (run in 1:num_runs) {
    cat(sprintf("Tissue %s: Run %d of stochastic search for MAX influence...\n", tissue, run))
    result_max <- stochastic_search("max", candidate_genes_tissue, graphs_tissue, lambda, epsilon, s_lambda)
    cat(sprintf("Tissue %s: Run %d of stochastic search for MIN influence...\n", tissue, run))
    result_min <- stochastic_search("min", candidate_genes_tissue, graphs_tissue, lambda, epsilon, s_lambda)
    
    # Create rows: Tissue, Metric, Gene Set (as concatenated string)
    max_row <- c(tissue, result_max$metric, paste(result_max$gene_set, collapse = ", "))
    min_row <- c(tissue, result_min$metric, paste(result_min$gene_set, collapse = ", "))
    
    # Append each run's row to the corresponding list.
    max_runs_list[[run]][[length(max_runs_list[[run]]) + 1]] <- max_row
    min_runs_list[[run]][[length(min_runs_list[[run]]) + 1]] <- min_row
    
    # Update best result for tissue.
    if (result_max$metric > best_max_metric) {
      best_max_metric <- result_max$metric
      best_max_row <- max_row
    }
    if (result_min$metric < best_min_metric) {
      best_min_metric <- result_min$metric
      best_min_row <- min_row
    }
  }
  
  # Store best rows for this tissue.
  max_best_list[[length(max_best_list) + 1]] <- best_max_row
  min_best_list[[length(min_best_list) + 1]] <- best_min_row
  
  cat(sprintf("Tissue %s processed. Best MAX metric: %s, Best MIN metric: %s\n", 
              tissue, best_max_metric, best_min_metric))
}

# -----------------------------
# Convert lists to data frames for output.
# Each run sheet: first column Tissue, second column Metric, third column Gene Set.
# -----------------------------
convert_to_df <- function(results_list) {
  # results_list is a list of lists (each inner list is a row vector)
  # First, determine the maximum number of columns (should be 3)
  max_ncols <- max(sapply(results_list, length))
  fill_row <- function(x, target_length) {
    length(x) <- target_length
    x
  }
  mat <- do.call(rbind, lapply(results_list, fill_row, target_length = max_ncols))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- c("Tissue", "AgeIndex", "Gene_Set")
  df
}

# For each run in max mode.
max_run_dfs <- list()
for (i in 1:num_runs) {
  max_run_dfs[[paste0("run", i)]] <- convert_to_df(max_runs_list[[i]])
}

# For each run in min mode.
min_run_dfs <- list()
for (i in 1:num_runs) {
  min_run_dfs[[paste0("run", i)]] <- convert_to_df(min_runs_list[[i]])
}

# Convert the best lists to data frames.
max_best_df <- convert_to_df(max_best_list)
min_best_df <- convert_to_df(min_best_list)

# Combine into a sheet list (adding a "best" sheet).
max_sheet_list <- c(max_run_dfs, list(best = max_best_df))
min_sheet_list <- c(min_run_dfs, list(best = min_best_df))

# -----------------------------
# Write the Excel output files.
# -----------------------------
output_folder <- file.path("..", "data", "age_indexes", "excel_files")
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

max_excel_file <- file.path(output_folder, "a1_max_age_indexes.xlsx")
min_excel_file <- file.path(output_folder, "a1_min_age_indexes.xlsx")

write_xlsx(max_sheet_list, path = max_excel_file)
write_xlsx(min_sheet_list, path = min_excel_file)

end_time_total <- Sys.time()
total_time_sec <- as.numeric(difftime(end_time_total, start_time_total, units = "secs"))
cat(sprintf("\nTotal runtime: %.3f seconds\n", total_time_sec))
