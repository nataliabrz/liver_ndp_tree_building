# Script for generation trees from NDP output

#############################################
## -- Receive command line args --
args <- commandArgs(trailingOnly = TRUE)
args

# Check for command line args, else set defaults
if (length(args) >= 3) {
  patient_id <- args[1]
  ndp_input_dir <- args[2]
  script_dir <- args[3]
} else {
  stop("Script needs at least 4 arguments.", call.=FALSE)
}

if (length(args) >= 4) {
  tree_out_dir <- args[4]
} else {
  tree_out_dir <- ndp_input_dir
}

if (length(args) >= 5) {
  file_snv <- args[5]
} else {
  snv_dir <- dirname(dirname(ndp_input_dir))
  file_snv <- paste(snv_dir,
                    paste0(patient_id, "_/bb_pass_snvs_all.csv"),
                    sep = "/")
  message("Using default snv file location of ", file_snv)
}

if (length(args) >= 6) {
  min_vaf_threshold <- args[6]
} else {
  min_vaf_threshold <- 0.10
  message("Using default min VAF threshold of 0.10")
}

if (length(args) >= 7) {
  min_mut_count <- args[7]
} else {
  min_mut_count <- 50
  message("Using default min mutation threshold of 50")
}

### debugging
if (FALSE) {
  script_dir <- "scripts/"
  patient_id <- "PD51606"
  ndp_input_dir <- "outputs/ndp_results/PD51606"
  snv_dir <- "data/filtered_calls/snv/"
  file_snv <- "data/filtered_calls/snv/PD51606_bb_pass_snvs_all.csv"
  tree_out_dir <- "outputs/trees/PD51606"
  min_vaf_threshold <- 0.10 #default 0.10
  min_mut_count <- 50 # default 50
}

# debug 20240617
print("Script dir-----")
print(script_dir)


### Load required Libraries
library(RColorBrewer)
library(philentropy)
library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(igraph)
library(ape)
library(stringr)
source(paste0(script_dir, "/utils/ndp_tree_functions.R"))
source(paste0(script_dir, "/utils/phylogeny_functions.R"))

### Get clean mut spectrum
mut_spec_by_clust <- extract.mut.spec(ndp_input_dir)

### Assign mut cluser ID
x_cluster_and_samples <- cluster.append(snvdir =  dirname(file_snv),
                                        ndpdir = ndp_input_dir,
                                        mut.spec.by.clust = mut_spec_by_clust,
                                        id.patient = patient_id)

### Filter to remove clusters not present in any sample above min vaf threshold
x_cluster_and_samples <- x_cluster_and_samples %>%
  filter(
    apply(
      x_cluster_and_samples[, 1:(ncol(x_cluster_and_samples) - 2)],
      1,
      function(x) any(x >= min_vaf_threshold)
    )
  )

save(x_cluster_and_samples,
     mut_spec_by_clust,
     file = paste(ndp_input_dir,
                  paste(patient_id, "x.clusters.RData", sep = "."),
                  sep = "/"))

### get number of mutations for branch lengths
mut_df <- extract.mutation.counts(x.cluster_and_samples = x_cluster_and_samples)
num_muts <- x_cluster_and_samples$no.of.mutations.assigned
names(num_muts) <- x_cluster_and_samples$cluster_id

### extract 95% confidence interval data
x_clusters_centiles <-
  extract.95.CI(x.cluster_and_samples = x_cluster_and_samples,
                ndpdir = ndp_input_dir)

### format median VAF tbl
x_cluster_and_samples_final <-
  create.med.vaf.tbl(x.cluster_and_samples = x_cluster_and_samples)

### asssign clusters to annotated mutations
snv_tbl <- fread(file_snv) %>%
  dplyr::mutate(pos_id = paste(Chrom, Pos, sep = "_"))
cluster_file <- paste(ndp_input_dir, "clust_assign_posthoc.csv", sep = "/")
cluster_tbl <- fread(cluster_file) %>%
  dplyr::mutate(pos_id = paste(chrom, pos, sep = "_")) %>%
  dplyr::select(-mut_id, -chrom, -pos)
assigned_tbl <- left_join(snv_tbl, cluster_tbl, by = "pos_id") %>%
  dplyr::rename("cluster_id" = clust_assign)
write.csv(assigned_tbl,
          file = paste(dirname(file_snv),
                       paste0(patient_id, "_ndp_assigned_muts.csv"),
                       sep = "/"), quote = FALSE, row.names = FALSE)

### save intermediate r data
save(mut_df,
     num_muts,
     x_clusters_centiles,
     x_cluster_and_samples_final,
     file = paste(ndp_input_dir,
                  paste(patient_id, "x.clusters.unique.fullrun.RData",
                        sep = "."),
                  sep = "/"))

### generate all combinations of cluster pairs that pass thresholds
get_cluster_pairs <- function(x_clusters_unique_real_merged,
                              min_num_mut = min_mut_count,
                              min_contrib = min_vaf_threshold) {
  target_clust <- rownames(x_clusters_unique_real_merged)

  heat_cut <- x_clusters_unique_real_merged[
    which(rownames(x_clusters_unique_real_merged) %in% target_clust),
  ] > min_contrib

  pair_list <- list()
  for (i in seq_len(ncol(heat_cut))) {
    inc_clust <- names(which(heat_cut[, i]))
    if (length(inc_clust) < 2) { next }
    clust_combs <- t(combn(inc_clust, m = 2))
    this_pairs <- tibble(cl1 = clust_combs[, 1], cl2 = clust_combs[, 2])
    pair_list[[i]] <- this_pairs
  }
  all_pairs <- do.call("rbind", pair_list) %>% distinct()
  all_pairs <- all_pairs %>%
    mutate(pair_id = sprintf("%s.%s", cl1, cl2),
           pair_rev_id = sprintf("%s.%s", cl2, cl1))
  all_pairs <- all_pairs %>% filter(!(pair_id %in% pair_rev_id))

  return(all_pairs %>% dplyr::select(-pair_id, -pair_rev_id))
}

### plot pigeonhold plots
centiles_tbl_all <- create.diag.plots(
  x.cluster_and_samples = x_cluster_and_samples_final,
  x.clusters.centiles = x_clusters_centiles,
  outdir = ndp_input_dir,
  id.patient = patient_id,
  num.muts = num_muts,
  min_num_mut = min_mut_count
)


clust_pairs <- unique(centiles_tbl_all[, c("cl1", "cl2")])

### determine tree branches
x_branches <- infer.tree.branches(clust_pairs = clust_pairs,
                                  centiles_tbl_all = centiles_tbl_all,
                                  min_contrib = min_vaf_threshold)

write.csv(x_branches,
          file = paste(tree_out_dir,
                       paste(patient_id, 'initial_branches.csv',sep='.'),
                       sep='/'),
          row.names = FALSE, quote = FALSE)


# decide which branches to keep
x_branches_keep <- keep(x_branches)
x_branches_keep <- x_branches_keep[!duplicated(x_branches_keep[,c('From','To')]),]
idx <- which(x_branches_keep$From %in% 'PD')
ids.clusters <- unique(sort(as.character(x_branches_keep$From)[-idx]))
idx <- which(ids.clusters %in% x_branches_keep$To[-idx])

if (length(idx) > 0) {
  ids.clusters <- ids.clusters[-idx]
}

x_branches_keep$Evidence <- NA
rownames(x_branches_keep) <- seq(1,dim(x_branches_keep)[1])

# assign strength of evidence of clonal relationships
x_branches <- as.data.frame(x_branches)
for (i in seq_len(nrow(x_branches_keep))) {
  idx <- which(
    x_branches$From %in% x_branches_keep$From[i] &
    x_branches$To %in% x_branches_keep$To[i]
  )
  x_branches_keep$Evidence[i] = as.character(x_branches$Evidence[idx])
  x_branches_keep$Assign = "auto"
}
save(x_branches_keep,
     centiles_tbl_all,
     file = paste(ndp_input_dir,
                 paste(patient_id, "x_branches_pre_tiebreak.RData", sep = "."),
                 sep = "/"))

### break ties if multiple potential paths, keeping higher VAF
x_branches_keep_final <- break.ties(x.branches.keep = x_branches_keep,
                                    centiles_tbl_all = centiles_tbl_all)

### save intermediate r data
save(mut_df, x_branches_keep_final,
     file = paste(ndp_input_dir,
                  paste(patient_id, "x_branches_keep.RData", sep = "."),
                  sep = "/"))

### prep tree building data
x_branches_keep_final$To <- gsub("Cl.",
                                 "",
                                 x_branches_keep_final$To,
                                 ignore.case = TRUE)
x_branches_keep_final$From <- gsub("Cl.",
                                   "",
                                   x_branches_keep_final$From,
                                   ignore.case = TRUE)
muts_per_cluster <- mut_df

write.csv(muts_per_cluster,
          file = paste(tree_out_dir,
                       paste(patient_id, "branch_lengths.csv", sep = "."),
                       sep = "/"),
          row.names = FALSE)

### create a new ape tree object
tr <- convert_edges_to_phylo(tree_tbl =  x_branches_keep_final,
                             muts_per_cluster =   muts_per_cluster)

### save treefile
write.tree(tr,
           file = paste(tree_out_dir,
                        paste(patient_id, "cluster_tree.phylo", sep = "."),
                        sep = "/"))

### plot tree
pdf(paste(tree_out_dir,
          paste(patient_id, "cluster_tree_rewrite.pdf", sep = "."),
          sep = "/"))

par(xpd = TRUE, mar = c(4.1, 2.1, 2.1, 6.1), family = "sans")
draw_nice_tree(tr, show_internal_nodes = TRUE,
               edge_strength = tr$edge_strength,
               edge_assign = tr$edge_assign,
               patient_id = substr(patient_id, 1, 7))
highlight_tree_nodes(tr, pch = 23, bg = "red", show_internal_nodes = TRUE)
dev.off()
