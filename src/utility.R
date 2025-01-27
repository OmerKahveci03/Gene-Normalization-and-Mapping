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