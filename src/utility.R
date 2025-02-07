# utility.R
# This script contains functions to be sourced by other scripts

influence_g1_g2 <- function(graph, g1, g2, lambda = 1) {
  d <- distances(graph, v = g1, to = g2)
  return(exp(-d * lambda))
}

influence_S_g <- function(graph, S, g, lambda = 1) {
  d_min <- min(distances(graph, v = S, to = g))
  return(exp(-d_min * lambda))
}

influence_S_V <- function(graph, S, lambda = 1) {
  all_genes  <- V(graph)$name
  inf_values <- sapply(all_genes, function(v) influence_S_g(graph, S, v, lambda))
  return(mean(inf_values))  # sum(...) / length(...)
}

# Returns influence_list to be used by age_based_influence
compute_influence_list <- function(age_groups,
                                   tissue_name,
                                   gene_set,
                                   lambda,
                                   data_dir = "..") {
  # Create an empty numeric vector to store influence values
  influence_list <- numeric(length(age_groups))
  
  for (k in seq_along(age_groups)) {
    ag <- age_groups[k]
    
    # Build the mapping_file path
    mapping_file <- file.path(data_dir, "data", "mapped",
                              as.character(ag),
                              paste0("mapped_", tissue_name, "_", ag, ".csv"))
    
    # Check if file exists
    if (!file.exists(mapping_file)) {
      cat(sprintf("Age group %d: File does not exist: %s\n", ag, mapping_file))
      influence_list[k] <- NA
      next
    }
    
    edges <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)
    if (nrow(edges) == 0) {
      cat(sprintf("Age group %d: No edges found for %s.\n", ag, tissue_name))
      influence_list[k] <- NA
      next
    }
    
    g <- graph_from_data_frame(edges, directed = FALSE)
    
    # Calculate Influence(S, V) for this age group
    result <- influence_S_V(g, gene_set, lambda)
    influence_list[k] <- result
  }
  
  return(influence_list)
}

# returns a value between -1 and 1 indicating age based influence
age_based_influence <- function(influence_list, epsilon = 0.005) {
  val <- 0
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
  
  val <- val / 15
  return(val)
}