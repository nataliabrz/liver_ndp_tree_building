library(data.table)
library(dplyr)

source("scripts/utils/pc8_Allocate_new_muts.R")

ndp_out_path <- "outputs/ndp_results/"
snv_path <- "data/filtered_calls/snv/"
indel_path <- "data/filtered_calls/indel/"

samples <- list.files(ndp_out_path, pattern = "^PD")

for (i in seq_along(samples)){

  snvs <- read.table(paste0(snv_path, samples[i], "_bb_pass_snvs_all.csv"),
                     sep = ",", header = TRUE)

  clust_assign <- read.table(paste0(ndp_out_path, "clust_assign_posthoc.csv"),
                             sep = ",", header = TRUE)

  snvs$pos_id <- paste0(snvs$Chrom, "_", snvs$Pos)
  clust_assign$pos_id <- paste0(clust_assign$chrom, "_", clust_assign$pos)

  if (length(unique(snvs$mut_id)) > nrow(clust_assign)) {
    missing <- setdiff(snvs$pos_id, clust_assign$pos_id)
    snvs <- snvs[-c(which(snvs$pos_id %in% missing)), ]
  }

  clust_assign$mut_no <- clust_assign$mut_id
  clust_assign$mut_id <- unique(snvs$mut_id)
  clust_assign <- clust_assign[, c("mut_id", "pos_id", "clust_assign")]
  colnames(clust_assign) <- c("mut_id", "pos_id_clust", "cluster_id")

  snvs_temp <- merge(snvs, clust_assign, by = "mut_id")


  if (sum(snvs_temp$pos_id != snvs_temp$pos_id_clust) == 0) {

    snvs_assigned <- snvs_temp %>%
      dplyr::select(-c(pos_id_clust))

  }else {
    print("position IDs don't match, check the cluster assignments")
  }

  rm(snvs_temp)

  ### now assign indels to clusters

  ls_merge_out <- readRDS(paste0(ndp_out_path, samples[i],
                                 "/object_lymph_merge_ls.dat"))
  post_ls_dist <- readRDS(paste0(ndp_out_path, samples[i],
                                 "/object_post_cluster_pos.dat"))
  new_y <- read.table(paste0(indel_path, samples[i], "_ndp_alt_bb_flt.csv"),
                      sep = ",", header = TRUE)
  new_n <- read.table(paste0(indel_path, samples[i], "_ndp_depth_bb_flt.csv"),
                      sep = ",", header = TRUE)

  new_y <- as.matrix(new_y[, 7:length(colnames(new_y))])
  new_n <- as.matrix(new_n[, 7:length(colnames(new_n))])

  new_cluster <- allocate.new.muts.to.cluster(ls_merge_out, post_ls_dist,
                                              new_y, new_n)

  context <- read.table(paste0(indel_path, samples[i],
                               "_mut_context_GRCh38.txt"),
                        sep = "\t", header = TRUE)

  context$clust_assign <- new_cluster
  context$mut_id <- rownames(context)
  cluster_file <- context[, c("clust_assign", "mut_id", "chrom", "pos")]

  file_indel <- paste0(indel_path, samples[i], "_bb_pass_indels_all.csv")

  ### assign clusters to annotated mutations
  indel_tbl <- fread(file_indel) %>%
    dplyr::mutate(pos_id = paste(Chrom, Pos, sep = "_"))
  cluster_tbl <- fread(cluster_file) %>%
    dplyr::mutate(pos_id = paste(chrom, pos, sep = "_")) %>%
    dplyr::select(-mut_id, -chrom, -pos)
  assigned_tbl <- left_join(indel_tbl, cluster_tbl, by = "pos_id") %>%
    dplyr::rename("cluster_id" = clust_assign)

  # combine snvs and indels into single dataframe
  assigned_all <- rbind(snvs_assigned, assigned_tbl)

  write.csv(assigned_all,
            file = paste0(ndp_out_path, samples[i], "/", samples[i],
                          "_ndp_assigned_snvs_and_indels.csv"),
            quote = FALSE, row.names = FALSE)

}