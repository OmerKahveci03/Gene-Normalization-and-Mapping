## How to set up
    - Download and setup RStudio
    - Change working directory to src
    - Download raw GTEX data from 'https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression' and untar it
    - See scripts below and run whichever you like. Start with 'normalize_all_tissues.R' and 'gene_networks.R'

## Directories
### data/raw
    -Store gtex data here
    -Naming convention: 'gene_tpm_{tissue name}.gctannotated.txt'

### data/normalized
    -Stores normalized data created by 'normalize_tissue.R' or 'normalize_all_tissues.R'
    -Naming convention: 'normalized_{tissue name}.csv'

### data/mapped
    - Contains six age group directories
    - Each age group directory contains a 'mapped_{tissue name}_{age group}.csv' file
    - Data created by 'gene_networks.R'

### src
    -Contains all R scripts

## R Scripts
### normalize_tissue.R
    -Creates a 'normlaized_{tissue name}.csv' file for a single input file
    -Removes the top epsilon % gene expression values
    -Normalizes values to be between 0 and 1
    -Takes the average of all expression values per age group
    -Creates normalized data as output

### normalize_all_tissues.R
    -Creates a 'normalized_{tissue name}.csv' file for all tissue data.
    -Same logic as 'normalize_tissue.R'

### gene_networks.R
    - Creates a 'mapped_{tissue name}_{age group}.csv' file for each age group for each tissue.
    - Uses normalized data as input
    - Reads the Trrust gene network and filters out genes with a normalized expression value less than threshhold value

### age_based_influence.R
    - Reads mapped data as input
    - Has a configurabe gene_set and tissue_name variable
    - Prints to the screen the age based influence of the gene set on the entire network, for the given tissue
    - Also creates a plot

### algorithm_0.R
    - Reads mapped data and trrust data as input
    - Has a configurable tissue_name variable
    - Iterates through all unique gene symbols in trrust column A
    - Finds the gene set within the given tissue name's gene network with the highest and lowest collective age-based influence value
    - Creates plots of the influence values per age group for the sets.