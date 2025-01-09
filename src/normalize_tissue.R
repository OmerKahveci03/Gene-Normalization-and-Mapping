# normalize_tissue.R
# This script will take in a .txt file of gene expression data.
# Everything is separated by tabs
# line three contains gene symbol names
# each line that comes after contains gene expression values for each subject
# The last four elements for each test subject is metadata.
# The second to last element is the subject's age group.
# The other three pieces of data will be ignored
# Age groups will go from 20-29, 30-39, to 70-79. Doesn't go higher or lower
# we will create a normalized_tissue_name.csv file in the output directory at the end

base_dir <- normalizePath("../")
input_dir <- file.path(base_dir, "data/raw")
output_dir <- file.path(base_dir, "data/normalized")

input_file <- "gene_tpm_brain_cortex.gctannotated.txt"
tissue_name <- sub("gene_tpm_(.*)\\.gctannotated\\.txt", "\\1", basename(input_file))

# configurable variable representing the amount of expression values we will ignore
# In this case it is the top 5%
epsilon <- 0.05


# To-do: read the raw data into a table or data frame

# To-do: Go through each column (gene) and create a list of the gene exps in ascending order
  # each gene exp in the list should have the subject's 'age group' stored as well.
  # delete the top epsilon% of the gene exp from the list
  # Then, take the min and max expression values in the list (post epsilon removal)
  # apply the formula v = (v-min)/(max-min) to each expression value (v)
  # Next, for each age group, take the average of all exp values. (there will be 6 age groups)
  # Add these values into a new table. The new table has 6 rows, one for each age group.
  # The columns for this new table will be genes.

# After we are finished, create the normalized file in the output directory.
  # Columns: Gene Symbols. Rows: The six age groups