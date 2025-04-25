#!/usr/bin/env Rscript

# null_algorithm.R
#  • build null distribution by sampling k random genes (k from a1_max best run)
#  • compute observed best‐run metrics from a1_max_age_indexes.xlsx and a1_min_age_indexes.xlsx
#  • compute Z scores: (v_obs – null_mean) / null_sd
#  • write two sheets: null_model and Statistical Significance

library(igraph)
library(readxl)
library(writexl)

# -----------------------------
# Paths & parameters
# -----------------------------
max_best_file <- file.path("..","data","age_indexes","excel_files",
                           "a1_max_age_indexes.xlsx")
min_best_file <- file.path("..","data","age_indexes","excel_files",
                           "a1_min_age_indexes.xlsx")
best_sheet    <- "best"

output_file   <- file.path("..","data","age_indexes","excel_files",
                           "null_age_indexes.xlsx")
n_runs        <- 50
lambda        <- 1
epsilon       <- 0.005

# -----------------------------
# Load TRRUST
# -----------------------------
trrust_file <- file.path("..","data","raw","trrust_rawdata.human.tsv")
if (!file.exists(trrust_file)) stop("TRRUST not found: ", trrust_file)
trrust       <- read.delim(trrust_file, sep="\t", header=TRUE,
                           stringsAsFactors=FALSE)
trrust_genes <- unique(c(trrust[[1]], trrust[[2]]))

# -----------------------------
# Helper functions
# -----------------------------
initialize_distances <- function(g) {
  if (is.null(g)) return(NULL)
  rep(Inf, vcount(g))
}
update_influence_for_graph <- function(g, d_cur, cand, lambda) {
  nodes <- V(g)$name
  if (!(cand %in% nodes)) {
    return(list(d_new=d_cur,
                influence_new=mean(exp(-lambda*d_cur))))
  }
  d_cand <- distances(g, v=cand, to=nodes)
  d_new  <- pmin(d_cur, as.vector(d_cand))
  list(d_new=d_new,
       influence_new=mean(exp(-lambda*d_new)))
}
age_based_influence <- function(inf_list, eps=0.005) {
  s <- 0
  for (i in 2:5) for (j in (i+1):6) {
    if (is.na(inf_list[i])||is.na(inf_list[j])) next
    if (inf_list[i]+eps < inf_list[j]) s <- s+1
    else if (inf_list[j]+eps < inf_list[i]) s <- s-1
  }
  s/15
}

# -----------------------------
# Read best‐run outputs (Metric column)
# -----------------------------
max_best_df <- read_excel(max_best_file, sheet=best_sheet)
min_best_df <- read_excel(min_best_file, sheet=best_sheet)

tissues   <- max_best_df$Tissue
ks        <- sapply(max_best_df$Gene_Set, function(s)
  length(strsplit(s, ",\\s*")[[1]]))
names(ks) <- tissues

v_max_obs <- max_best_df$Metric
v_min_obs <- min_best_df$Metric

# -----------------------------
# Prepare null‐model storage
# -----------------------------
results <- data.frame(
  Tissue = tissues,
  Mean   = NA_real_,
  SD     = NA_real_,
  matrix(NA_real_, nrow=length(tissues), ncol=n_runs,
         dimnames=list(NULL, paste0("run",1:n_runs))),
  stringsAsFactors=FALSE
)

# -----------------------------
# Build null distribution
# -----------------------------
for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  k      <- ks[tissue]
  cat(sprintf("=== %s: sampling %d genes/run ===\n", tissue, k))
  
  # load graphs & mapped_src
  graphs    <- list()
  mapped_src <- character(0)
  for (ag in c(20,30,40,50,60,70)) {
    fn <- file.path("..","data","mapped",as.character(ag),
                    sprintf("mapped_%s_%d.csv", tissue, ag))
    if (!file.exists(fn)) { graphs[[as.character(ag)]]<-NULL; next }
    edges <- read.csv(fn, stringsAsFactors=FALSE)
    if (nrow(edges)==0) { graphs[[as.character(ag)]]<-NULL; next }
    g <- graph_from_data_frame(edges, directed=TRUE)
    graphs[[as.character(ag)]] <- g
    mapped_src <- union(mapped_src, unique(edges[[1]]))
  }
  
  # candidate genes from TRRUST ∩ mapped_src
  cands <- intersect(trrust_genes, mapped_src)
  if (length(cands) < k) {
    warning(sprintf("%s: only %d candidates (<k=%d)—filling NA", 
                    tissue, length(cands), k))
    results[i, paste0("run",1:n_runs)] <- NA_real_
    next
  }
  
  # 50 null runs
  metrics <- numeric(n_runs)
  for (r in 1:n_runs) {
    sel <- sample(cands, k)
    ds  <- lapply(graphs, initialize_distances)
    for (gname in sel) {
      for (j in seq_along(graphs)) {
        gr <- graphs[[j]]
        if (!is.null(gr)) {
          res <- update_influence_for_graph(gr, ds[[j]], gname, lambda)
          ds[[j]] <- res$d_new
        }
      }
    }
    infs      <- mapply(function(g,d) {
      if (is.null(g)) return(NA)
      mean(exp(-lambda*d))
    }, graphs, ds)
    metrics[r] <- age_based_influence(infs, epsilon)
  }
  
  # store
  results[i, paste0("run",1:n_runs)] <- metrics
  results$Mean[i] <- mean(metrics, na.rm=TRUE)
  results$SD[i]   <- sd(metrics, na.rm=TRUE)
}

# -----------------------------
# Statistical Significance (Z‐scores)
# -----------------------------
obs_df <- data.frame(
  Tissue  = tissues,
  Max_Obs = v_max_obs,
  Min_Obs = v_min_obs,
  stringsAsFactors=FALSE
)
sig_df <- merge(results, obs_df, by="Tissue", all.x=TRUE)
sig_df$Max_Obs <- as.numeric(sig_df$Max_Obs)
sig_df$Min_Obs <- as.numeric(sig_df$Min_Obs)
sig_df$Mean    <- as.numeric(sig_df$Mean)
sig_df$SD      <- as.numeric(sig_df$SD)

sig_df$`Max Data Z Score` <- (sig_df$Max_Obs - sig_df$Mean) / sig_df$SD
sig_df$`Min Data Z Score` <- (sig_df$Min_Obs - sig_df$Mean) / sig_df$SD

stat_sig <- sig_df[, c("Tissue", "Min Data Z Score", "Max Data Z Score")]

# -----------------------------
# Write both sheets
# -----------------------------
write_xlsx(
  list(
    null_model               = results,
    `Statistical Significance` = stat_sig
  ),
  path=output_file
)

cat("Written null_model + Statistical Significance to", output_file, "\n")
