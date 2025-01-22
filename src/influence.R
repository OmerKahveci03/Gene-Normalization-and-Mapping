# influence.R
# Given tissue name and set of genes
# Find the age-based-influence value

# Create function influence(g1, g2) where g1 and g2 are genes
# influence(g1, g2) = e^(-dist(g1, g2) * lambda)
# dist(g1, g2) is the shortest path length between g1 and 

# influence(S, g) = e^(-min(dist(gi, g)) * lambda)
# S is a set of genes, gi is a gene in the set.
# min() means find the gene in the set with the shortest distance

# influence(S, V) = 1/(size of V) * (sum of all inf(S, vi))
# V is the set of ALL genes, vi is the individual gene in V
# (size of V) is the total number of genes

# Connectivity of genes is found in the mapping_dir.
  # in mapping_dir, there are six directories named {age group}
    # {age group} is 20, 30, 40, 50, 60, or 70
    # each {age group} directory has a mapped csv file for all tissues
    # Naming convention is mapped_{tissue name}_{age group}.csv


base_dir    <- normalizePath("../src")
mapping_dir <- file.path(base_dir, "../data/raw")

tissue_name <- "brain_cortex"
lambda <- 1
# gene_set <- {'AATF', 'ABL1', 'AES'}

# Let's start small. Give me the influence(S, V) for the tissue brain_cortex at age group 20 and the commented out gene_set
# influence(S, V) will be done by running influence(S, g) (for g in V)
# influence(S,g) will be done by running influence(g1, g2) (for g1 in S, and g2 is g)
# print the result to the screen using cat.