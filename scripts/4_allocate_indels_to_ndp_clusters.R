
library(data.table)
library(dplyr)

source("scripts/utils/pc8_Allocate_new_muts.R")

ndp_out_path <- "outputs/ndp_results/"
snv_path <- "data/filtered_calls/snv/"
indel_path <- "data/filtered_calls/indel/"

samples <- list.files(ndp_out_path, pattern = "^PD")

for(i in seq_along(samples)){

ls.merge.out <- readRDS(paste0(ndp_out_path,samples[i],"/object_lymph.merge.ls.dat"))
post.ls.dist <- readRDS(paste0(ndp_out_path,samples[i],"/object_post.cluster.pos.dat"))
new.y <- read.table(paste0(indel_path,samples[i],"_ndp_alt_bb_flt.csv"),sep=",",header=T)
new.N <- read.table(paste0(indel_path,samples[i],"_ndp_depth_bb_flt.csv"),sep=",",header=T)

new.y <- as.matrix(new.y[,7:length(colnames(new.y))])
new.N <- as.matrix(new.N[,7:length(colnames(new.N))])

new.cluster <- allocate.new.muts.to.cluster(ls.merge.out, post.ls.dist, new.y, new.N)

context <- read.table(paste0(indel_path,samples[i],"_mut_context_GRCh38.txt"),sep="\t",header=T)
context$clust_assign <- new.cluster
context$mut_id <- rownames(context)
cluster.file <- context[,c("clust_assign","mut_id","chrom","pos")]

file.indel <- paste0(indel_path,samples[i],"_bb_pass_indels_all.csv")

### assign clusters to annotated mutations
indel_tbl <- fread(file.indel) %>%
  dplyr::mutate(pos_id = paste(Chrom, Pos, sep = "_"))
cluster_tbl <- fread(cluster.file) %>%
  dplyr::mutate(pos_id = paste(chrom, pos, sep = "_")) %>%
  dplyr::select(-mut_id, -chrom, -pos)
assigned_tbl <- left_join(indel_tbl, cluster_tbl, by = "pos_id") %>%
  dplyr::rename("cluster_id" = clust_assign)
write.csv(assigned_tbl, file = paste0(ndp_out_path,samples[i],"/",samples[i],"_ndp_assigned_indels.csv"), quote = F, row.names = F)

# combine snvs and indels into single dataframe
assigned_snvs <- read.table(paste0(ndp_out_path,samples[i],"/",samples[i],"_ndp_assigned_snvs.csv"),sep=",",header=T)
assigned_all <- rbind(assigned_snvs, assigned_tbl)
write.csv(assigned_all,file = paste0(ndp_out_path,samples[i],"/",samples[i],"_ndp_assigned_snvs_and_indels.csv"), quote = F, row.names = F)

}