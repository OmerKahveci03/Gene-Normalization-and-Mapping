## Directories
### data/raw
    -Store gtex data here
    -Naming convention: `gene_tpm_{tissue name}.gctannotated.txt`

### data/normalized
    -Stores normalized data created by `normalize_tissue.R` or `normalize_all_tissues.R`
    -Naming convention: `normalized_{tissue name}.csv`

### src
    -Contains all R scripts

## R Scripts
### normalize_tissue.R
    -Creates a normlaized_{tissue name}.csv file for a single input file
    -Removes the top epsilon % gene expression values
    -Normalizes values to be between 0 and 1
    -Takes the average of all expression values per age group
    -Creates normalized data as output

### normalize_all_tissues.R
    -Creates a normalized_{tissue name}.csv file for all tissue data.
    -Same logic as `normalize_tissue.R`