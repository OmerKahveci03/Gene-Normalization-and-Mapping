# age_based_influence.R

# 1) influence(g1, g2)   = e^(-dist(g1, g2) * lambda)
# 2) influence(S, g)     = e^(- (min_{gi in S} dist(gi, g)) * lambda)
# 3) influence(S, V)     = (1 / |V|) * sum(influence(S, v)) for v in V

# Given a tissue name and gene set, find the age-based influence and plot the results.

library(igraph)

# Configurable Values
tissue_name  <- "brain_cortex"
lambda       <- 1
gene_set     <- c("AATF", "ABL1", "AES")
epsilon      <- 0.005

# Don't Configure These
age_groups          <- c(20, 30, 40, 50, 60, 70)
age_based_influence <- 0
influence_list      <- numeric(length(age_groups))

# --------------------------
# Define the influence functions
# --------------------------
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

# --------------------------
# Loop over each age group and compute influence
# --------------------------
for (k in seq_along(age_groups)) {
  ag <- age_groups[k]
  mapping_file <- file.path("..", "data", "mapped",
                            as.character(ag),
                            paste0("mapped_", tissue_name, "_", ag, ".csv"))
  
  if (!file.exists(mapping_file)) {
    cat(sprintf("Age group %d: File does not exist: %s\n", ag, mapping_file))
    influence_list[k] <- NA
    next
  }
  
  edges <- read.csv(mapping_file, header = TRUE, stringsAsFactors = FALSE)
  # If no edges, skip or handle gracefully
  if (nrow(edges) == 0) {
    cat(sprintf("Age group %d: No edges found for %s.\n", ag, tissue_name))
    influence_list[k] <- NA
    next
  }
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  result            <- influence_S_V(g, gene_set, lambda)
  influence_list[k] <- result
  
  cat(sprintf("Age group %2d => Influence(S, V): %f\n", ag, result))
}

# --------------------------
# Compute age-based influence
# --------------------------
for (i in 2:5) {
  for (j in (i+1):6) {
    if (is.na(influence_list[i]) || is.na(influence_list[j])) {
      next
    }
    if (influence_list[i] + epsilon < influence_list[j]) {
      age_based_influence <- age_based_influence + 1
    } else if (influence_list[j] + epsilon < influence_list[i]) {
      age_based_influence <- age_based_influence - 1
    }
  }
}
age_based_influence <- age_based_influence / 15

# --------------------------
# Plot the influence_list
# --------------------------

# 1) Determine the y-axis upper bound (round up to nearest tenth).
max_val  <- max(influence_list, na.rm = TRUE)
# If all values are NA, max_val might be -Inf, so let's handle that:
if (is.infinite(max_val)) {
  # No valid data, set a default
  max_val <- 1e-3
}

# Round up to the next tenth: e.g. 0.06 -> 0.1, 0.24 -> 0.3
ymax_rounded <- ceiling(max_val * 10) / 10
# Ensure at least something > 0
if (ymax_rounded <= 0) {
  ymax_rounded <- 0.1
}

plot(
  x    = age_groups,
  y    = influence_list,
  type = "b",         # "b" => draws both points and lines
  pch  = 19,          # Solid circle points
  col  = "blue",
  xlab = "Age Group",
  ylab = "Influence(S, V)",
  main = "Influence Values by Age Group",
  ylim = c(0, ymax_rounded)  # from 0 to the "rounded up" maximum
)

cat(sprintf("\nFinal age_based_influence: %f\n", age_based_influence))