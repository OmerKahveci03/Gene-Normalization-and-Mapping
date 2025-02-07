# age_based_influence.R

# Libraries & Utilities
library(igraph)
source("utility.R")  # <-- Brings in compute_influence_list() and age_based_influence()

# Configurable Values
tissue_name  <- "brain_cortex"
lambda       <- 1
gene_set     <- c("AATF", "ABL1", "AES")
epsilon      <- 0.005

age_groups <- c(20, 30, 40, 50, 60, 70)

#Compute influence_list
influence_list <- compute_influence_list(
  age_groups   = age_groups,
  tissue_name  = tissue_name,
  gene_set     = gene_set,
  lambda       = lambda,
  data_dir     = ".."  # adjust this if needed
)

for (k in seq_along(age_groups)) {
  cat(sprintf("Age group %2d => Influence(S, V): %f\n", age_groups[k], influence_list[k]))
}

#Compute age-based influence
age_based_influence_val <- age_based_influence(influence_list, epsilon)

# Plot the influence_list
max_val <- max(influence_list, na.rm = TRUE)
if (is.infinite(max_val)) {
  max_val <- 1e-3
}

# Round up to the next tenth
ymax_rounded <- ceiling(max_val * 10) / 10
if (ymax_rounded <= 0) {
  ymax_rounded <- 0.1
}

plot(
  x    = age_groups,
  y    = influence_list,
  type = "b",
  pch  = 19,
  col  = "blue",
  xlab = "Age Group",
  ylab = "Influence(S, V)",
  main = "Influence Values by Age Group",
  ylim = c(0, ymax_rounded)
)

cat(sprintf("\nFinal age_based_influence_val: %f\n", age_based_influence_val))
