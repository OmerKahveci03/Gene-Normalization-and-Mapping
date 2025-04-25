# Total runtime: 14503.127 seconds
library(igraph)
library(writexl)

# -----------------------------
# Determine Tissue List from mapped files for age 20
# -----------------------------
tissue_dir <- file.path("..", "data", "mapped", "20")
mapped_files <- list.files(path = tissue_dir,
                           pattern = "^mapped_.*_20\\.csv$",
                           full.names = FALSE)
tissues <- unique(gsub("^mapped_(.*)_20\\.csv$", "\\1", mapped_files))
cat(sprintf("Found %d tissues: %s\n",
            length(tissues), paste(tissues, collapse = ", ")))

# -----------------------------
# Configurable Values
# -----------------------------
lambda     <- 1
epsilon    <- 0.005
age_groups <- c(20, 30, 40, 50, 60, 70)
s_lambda   <- 5   # Constant to scale improvement scores for stochastic
num_runs   <- 5   # Number of stochastic runs per tissue

# -----------------------------
# Start overall timer
# -----------------------------
start_time_total <- Sys.time()

# -----------------------------
# Load TRRUST file to get all unique gene symbols
# -----------------------------
trrust_file <- file.path("..", "data", "raw", "trrust_rawdata.human.tsv")
if (!file.exists(trrust_file)) stop("TRRUST file not found at: ", trrust_file)
trrust <- read.delim(trrust_file,
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     sep = "\t")
trrust_genes <- unique(c(trrust[[1]], trrust[[2]]))

# -----------------------------
# Core helper functions
# -----------------------------
initialize_distances <- function(g) {
  if (is.null(g)) return(NULL)
  rep(Inf, length(V(g)$name))
}

update_influence_for_graph <- function(g, d_current, candidate, lambda) {
  all_nodes <- V(g)$name
  if (!(candidate %in% all_nodes)) {
    return(list(
      d_new = d_current,
      influence_new = mean(exp(-lambda * d_current))
    ))
  }
  d_candidate <- as.vector(distances(g,
                                     v = candidate,
                                     to = all_nodes))
  d_new <- pmin(d_current, d_candidate)
  influence_new <- mean(exp(-lambda * d_new))
  list(d_new = d_new, influence_new = influence_new)
}

age_based_influence <- function(inf_list, epsilon = 0.005) {
  val <- 0
  for (i in 2:5) {
    for (j in (i + 1):6) {
      if (is.na(inf_list[i]) || is.na(inf_list[j])) next
      if (inf_list[i] + epsilon < inf_list[j]) {
        val <- val + 1
      } else if (inf_list[j] + epsilon < inf_list[i]) {
        val <- val - 1
      }
    }
  }
  val / 15
}

stochastic_search <- function(mode = c("max", "min"),
                              candidate_genes, graphs,
                              lambda, epsilon, s_lambda) {
  mode <- match.arg(mode)
  gene_set <- character(0)
  current_d <- lapply(graphs, initialize_distances)
  current_inf <- mapply(function(g, d) {
    if (is.null(g)) return(NA)
    mean(exp(-lambda * d))
  }, graphs, current_d)
  current_metric <- age_based_influence(current_inf, epsilon)
  improved <- TRUE; iteration <- 0
  
  while (improved) {
    iteration <- iteration + 1
    cat(sprintf("\n%s Iteration %d\n",
                toupper(mode), iteration))
    remaining <- setdiff(candidate_genes, gene_set)
    
    candidate_list <- list()
    candidate_updates_list <- list()
    candidate_weights <- numeric()
    
    for (cand in remaining) {
      new_inf <- numeric(length(graphs))
      updates <- vector("list", length(graphs))
      
      for (i in seq_along(graphs)) {
        g <- graphs[[i]]
        if (is.null(g)) {
          new_inf[i] <- NA
          updates[[i]] <- NULL
        } else {
          res <- update_influence_for_graph(
            g, current_d[[i]], cand, lambda
          )
          updates[[i]] <- res$d_new
          new_inf[i] <- res$influence_new
        }
      }
      new_metric <- age_based_influence(new_inf, epsilon)
      
      good <- (mode == "max" && new_metric > current_metric) ||
        (mode == "min" && new_metric < current_metric)
      if (good) {
        if (current_metric == 0 && mode == "max") {
          imp_factor <- new_metric
          weight <- imp_factor * s_lambda
        } else {
          imp_factor <- if (mode == "max")
            new_metric / current_metric
          else
            current_metric / new_metric
          weight <- imp_factor^s_lambda
        }
        if (is.finite(weight) && weight > 0) {
          candidate_list <- c(candidate_list, cand)
          candidate_updates_list[[length(candidate_list)]] <- updates
          candidate_weights <- c(candidate_weights, weight)
        }
      }
    }
    
    # FIXED CONDITIONAL:
    # stop if no candidates or no positive weight
    if (length(candidate_list) == 0 || sum(candidate_weights) <= 0) {
      cat(sprintf("No candidate gene improves the metric further (%s).\n",
                  mode))
      improved <- FALSE
      break
    }
    
    sel <- sample(seq_along(candidate_list),
                  size = 1,
                  prob = candidate_weights)
    best_cand <- candidate_list[[sel]]
    best_updates <- candidate_updates_list[[sel]]
    
    # apply update
    gene_set <- c(gene_set, best_cand)
    for (i in seq_along(graphs)) {
      if (!is.null(graphs[[i]])) {
        current_d[[i]] <- best_updates[[i]]
      }
    }
    current_inf <- mapply(function(g, d) {
      if (is.null(g)) return(NA)
      mean(exp(-lambda * d))
    }, graphs, current_d)
    current_metric <- age_based_influence(current_inf, epsilon)
    
    cat(sprintf("  Chosen Gene: %s  New Metric: %f\n",
                best_cand, current_metric))
  }
  
  final_inf <- mapply(function(g, d) {
    if (is.null(g)) return(NA)
    mean(exp(-lambda * d))
  }, graphs, current_d)
  
  list(
    gene_set   = gene_set,
    metric     = current_metric,
    influences = final_inf
  )
}

# -----------------------------
# Add the brute-force greedy_search (algorithm 0)
# -----------------------------
greedy_search <- function(mode = c("max", "min"),
                          candidate_genes, graphs,
                          lambda, epsilon) {
  mode <- match.arg(mode)
  gene_set <- character(0)
  current_d <- lapply(graphs, initialize_distances)
  current_inf <- mapply(function(g, d) {
    if (is.null(g)) return(NA)
    mean(exp(-lambda * d))
  }, graphs, current_d)
  current_metric <- age_based_influence(current_inf, epsilon)
  
  improved <- TRUE; iteration <- 0
  while (improved) {
    iteration <- iteration + 1
    cat(sprintf("\n%s Iteration %d\n",
                toupper(mode), iteration))
    remaining <- setdiff(candidate_genes, gene_set)
    
    best_cand <- NA; best_metric <- current_metric
    best_updates <- NULL
    
    for (cand in remaining) {
      new_inf <- numeric(length(graphs))
      updates <- vector("list", length(graphs))
      
      for (i in seq_along(graphs)) {
        g <- graphs[[i]]
        if (is.null(g)) {
          new_inf[i] <- NA
          updates[[i]] <- NULL
        } else {
          res <- update_influence_for_graph(
            g, current_d[[i]], cand, lambda
          )
          updates[[i]] <- res$d_new
          new_inf[i] <- res$influence_new
        }
      }
      nm <- age_based_influence(new_inf, epsilon)
      
      if ((mode == "max" && nm > best_metric) ||
          (mode == "min" && nm < best_metric)) {
        best_metric <- nm
        best_cand   <- cand
        best_updates <- updates
      }
    }
    
    if (!is.na(best_cand) &&
        ((mode == "max" && best_metric > current_metric) ||
         (mode == "min" && best_metric < current_metric))) {
      
      gene_set      <- c(gene_set, best_cand)
      current_metric <- best_metric
      for (i in seq_along(graphs)) {
        if (!is.null(graphs[[i]])) {
          current_d[[i]] <- best_updates[[i]]
        }
      }
      cat(sprintf("  Chosen Gene: %s  New Metric: %f\n",
                  best_cand, current_metric))
    } else {
      cat(sprintf("No candidate gene improves the metric further (%s).\n",
                  mode))
      improved <- FALSE
    }
  }
  
  final_inf <- mapply(function(g, d) {
    if (is.null(g)) return(NA)
    mean(exp(-lambda * d))
  }, graphs, current_d)
  
  list(
    gene_set   = gene_set,
    metric     = current_metric,
    influences = final_inf
  )
}

# -----------------------------
# Prepare lists to store Excel results for each run (0…5) and best
# -----------------------------
total_runs    <- num_runs + 1
max_runs_list <- vector("list", total_runs)
min_runs_list <- vector("list", total_runs)
for (i in seq_len(total_runs)) {
  max_runs_list[[i]] <- list()
  min_runs_list[[i]] <- list()
}
max_best_list <- list()
min_best_list <- list()

# -----------------------------
# Loop over all tissues
# -----------------------------
for (tissue in tissues) {
  cat(sprintf("\n=============================\nProcessing Tissue: %s\n",
              tissue))
  
  # load all age‑group graphs for this tissue…
  graphs_tissue <- list()
  mapped_source_genes <- character(0)
  for (ag in age_groups) {
    f <- file.path("..", "data", "mapped", as.character(ag),
                   paste0("mapped_", tissue, "_", ag, ".csv"))
    if (!file.exists(f)) {
      graphs_tissue[[as.character(ag)]] <- NULL; next
    }
    edges <- read.csv(f, header=TRUE, stringsAsFactors=FALSE)
    if (nrow(edges) == 0) {
      graphs_tissue[[as.character(ag)]] <- NULL; next
    }
    g <- graph_from_data_frame(edges, directed=TRUE)
    graphs_tissue[[as.character(ag)]] <- g
    mapped_source_genes <- union(mapped_source_genes,
                                 unique(edges[[1]]))
  }
  
  candidate_genes <- intersect(trrust_genes, mapped_source_genes)
  if (length(candidate_genes) == 0) {
    cat(sprintf("Tissue %s: No candidates. Skipping.\n", tissue))
    next
  }
  
  # initialize best‑trackers
  best_max_metric <- -Inf
  best_min_metric <-  Inf
  best_max_row    <- NULL
  best_min_row    <- NULL
  
  # ----- Run 0: Brute‑force greedy_search -----
  cat(sprintf("Tissue %s: Run 0 (greedy) MAX...\n", tissue))
  r0_max <- greedy_search("max", candidate_genes,
                          graphs_tissue, lambda, epsilon)
  cat(sprintf("Tissue %s: Run 0 (greedy) MIN...\n", tissue))
  r0_min <- greedy_search("min", candidate_genes,
                          graphs_tissue, lambda, epsilon)
  
  row0_max <- c(tissue,
                r0_max$metric,
                paste(r0_max$gene_set, collapse = ", "))
  row0_min <- c(tissue,
                r0_min$metric,
                paste(r0_min$gene_set, collapse = ", "))
  
  max_runs_list[[1]][[length(max_runs_list[[1]]) + 1]] <- row0_max
  min_runs_list[[1]][[length(min_runs_list[[1]]) + 1]] <- row0_min
  
  if (r0_max$metric > best_max_metric) {
    best_max_metric <- r0_max$metric; best_max_row <- row0_max
  }
  if (r0_min$metric < best_min_metric) {
    best_min_metric <- r0_min$metric; best_min_row <- row0_min
  }
  
  # ----- Runs 1…5: Stochastic -----
  for (run in seq_len(num_runs)) {
    cat(sprintf("Tissue %s: Run %d (stochastic) MAX...\n",
                tissue, run))
    res_max <- stochastic_search("max", candidate_genes,
                                 graphs_tissue,
                                 lambda, epsilon, s_lambda)
    cat(sprintf("Tissue %s: Run %d (stochastic) MIN...\n",
                tissue, run))
    res_min <- stochastic_search("min", candidate_genes,
                                 graphs_tissue,
                                 lambda, epsilon, s_lambda)
    
    row_max <- c(tissue,
                 res_max$metric,
                 paste(res_max$gene_set, collapse = ", "))
    row_min <- c(tissue,
                 res_min$metric,
                 paste(res_min$gene_set, collapse = ", "))
    
    idx <- run + 1
    max_runs_list[[idx]][[length(max_runs_list[[idx]]) + 1]] <- row_max
    min_runs_list[[idx]][[length(min_runs_list[[idx]]) + 1]] <- row_min
    
    if (res_max$metric > best_max_metric) {
      best_max_metric <- res_max$metric; best_max_row <- row_max
    }
    if (res_min$metric < best_min_metric) {
      best_min_metric <- res_min$metric; best_min_row <- row_min
    }
  }
  
  # save best for this tissue
  max_best_list[[length(max_best_list) + 1]] <- best_max_row
  min_best_list[[length(min_best_list) + 1]] <- best_min_row
  
  cat(sprintf("Tissue %s done. Best MAX=%f, Best MIN=%f\n",
              tissue,
              best_max_metric,
              best_min_metric))
}

# -----------------------------
# Convert to data.frames and write Excel
# -----------------------------
convert_to_df <- function(rows) {
  max_n <- max(sapply(rows, length))
  mat <- do.call(rbind, lapply(rows, function(r) {
    length(r) <- max_n
    r
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- c("Tissue", "Metric", "Gene_Set")
  df
}

max_sheets <- list()
min_sheets <- list()
for (i in 0:num_runs) {
  max_sheets[[paste0("run", i)]] <-
    convert_to_df(max_runs_list[[i + 1]])
  min_sheets[[paste0("run", i)]] <-
    convert_to_df(min_runs_list[[i + 1]])
}
max_sheets[["best"]] <- convert_to_df(max_best_list)
min_sheets[["best"]] <- convert_to_df(min_best_list)

out_dir <- file.path("..", "data", "age_indexes", "excel_files")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_xlsx(max_sheets,
           path = file.path(out_dir, "a1_max_age_indexes.xlsx"))
write_xlsx(min_sheets,
           path = file.path(out_dir, "a1_min_age_indexes.xlsx"))

end_time_total <- Sys.time()
cat(sprintf("\nTotal runtime: %.3f seconds\n",
            as.numeric(difftime(end_time_total,
                                start_time_total,
                                units = "secs"))))
