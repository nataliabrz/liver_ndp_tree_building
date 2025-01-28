# Run Dirichlet clustering
# --
# Runs Dirichlet clustering and label switching on low input mutation data
#
# Args (command line):
#   1: Path to git repository
#      (ie. parent directory where repository was cloned into)
#   2: Path to CSV containing counts of alt base
#      (rows are variants, any columns starting with 'PD' are samples)
#   3: Path to CSV containing depth
#      (rows are variants, any columns starting with 'PD' are samples)
#   4: Path to table with mutation contexts
#      (generated with "1_pull_context.R" script)
#   5: Directory to output results
#   6: Number of burnin cycles, eg 10000
#
# If command line args are not provided, default file locations are used.
#
# --
# /// Author --- PETER CAMPBELL
# /// Edits --- Simon Brunner
# /// Edit date --- 11-MAY-2018
# /// Edits --- Stanley Ng
# /// Edit date --- 25-MAY-2019
# /// Edits --- Natalia Brzozowska
# /// Edit date --- 13-FEB-2023
#############################################
## -- Receive command line args --
args <- commandArgs(trailingOnly = TRUE)
args

# Check for command line args, else set defaults
if (length(args) == 0) {
  repo_location <- "scripts/utils/"
  fpath_alt_csv <- "data/filtered_calls/snv/PD51606_ndp_alt_bb_flt.csv"
  fpath_depth_csv <- "data/filtered_calls/snv/PD51606_ndp_depth_bb_flt.csv"
  fpath_mutant_locations <- file.path(
    "data",
    "filtered_calls",
    "snv",
    "PD51606_mut_context_GRCh38.txt"
  )
  sex_file <- "data/metadata/sex_file.txt"
  output_dir <- "outputs/ndp_results/PD51606"
  num_burnin <- 15000
} else if (length(args) == 7) {
  repo_location <- args[1]
  fpath_alt_csv <- args[2]
  fpath_depth_csv <- args[3]
  fpath_mutant_locations <- args[4]
  sex_file <- args[5]
  output_dir <- args[6]
  num_burnin <- as.numeric(args[7])
} else {
  stop("Script takes either exactly 0 or 7 arguments.", call. = FALSE)
}

# Report arguments
message(sprintf("Reading alt counts from %s", fpath_alt_csv))
message(sprintf("Reading depth from %s", fpath_depth_csv))
message(sprintf("Reading mutation contexts from %s", fpath_mutant_locations))
message(sprintf("Writing output to %s", output_dir))
message(sprintf("Number of burnin cycles is %s", num_burnin))

#############################################
## -- Library loading and sourcing --

library(RColorBrewer)
library(label.switching)
library(philentropy)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
library(here)
library(farver)


source(file.path(repo_location,
                 "Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R"))
source(file.path(repo_location,
                 "Mutation_Cluster_Dirichlet_posthoc_merge_fx.R"))

dir.create(output_dir, showWarnings = FALSE)

# Define params
number_iterations <- 25000
number_burn_in <- num_burnin

set.seed(28)

## define gender
f_list <- c("F", "f", "FEMALE", "female", "Female")
m_list <- c("M", "m", "MALE", "male", "Male")

#############################################
## -- Run clustering --

lymph_var <- read.table(fpath_alt_csv, sep = ",", header = TRUE,
                        stringsAsFactors = FALSE)

lymph_depth <- read.table(fpath_depth_csv, sep = ",", header = TRUE,
                          stringsAsFactors = FALSE)

mut_spec <- read.table(fpath_mutant_locations, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)

######### added by nb15:
#for males, double the coverage on X and Y
#for females, remove chr Y muts
sex_tbl <- fread(sex_file)
pd_columns <- colnames(lymph_depth)[grep("PD", colnames(lymph_depth))]
this_patient <- substr(pd_columns[1], 1, 7)
pat_sex <- sex_tbl[sex_tbl$patientID == this_patient]$sex

if (pat_sex %in% m_list) {
  cat("patient is male, doubling coverage on X and Y chromosomes")
  chrom_filter <- lymph_depth$chrom %in% c("chrX", "chrY")
  pd_columns <- grep("PD", colnames(lymph_depth))
  lymph_depth[chrom_filter, pd_columns] <-
    lymph_depth[chrom_filter, pd_columns] * 2
}else if (pat_sex %in% f_list) {
  cat("patient is female, removing muts on Y chromosome")
  lymph_var <- lymph_var[lymph_var$chrom != "chrY", ] # Remove Y chrom muts
  lymph_depth <- lymph_depth[lymph_depth$chrom != "chrY", ]
  mut_spec <- mut_spec[mut_spec$chrom != "chrY", ]
}

#########

lymph_y <- as.matrix(lymph_var[, which(grepl("PD", names(lymph_var)))])
lymph_n <- as.matrix(lymph_depth[, which(grepl("PD", names(lymph_depth)))])
lymph_gs <- muts.cluster.dirichlet.gibbs(C = 100, y = lymph_y, N = lymph_n,
                                         iter = number_iterations,
                                         iter.chunk = 100, outdir = output_dir)

short_names <- names(lymph_var)[which(grepl("PD", names(lymph_var)))]
short_names[2:length(short_names)] <-
  substring(short_names[2:length(short_names)], 8)

save.image(file.path(output_dir, "Rsession.dat"))

#############################################

# Generate convergence plots

pdf(file.path(output_dir, "Convergence_plots.pdf"),
    paper = "a4", width = 8, height = 10)

for (i in c(10, 50, 100, 200, 500, seq(1000, number_iterations, 1000))) {
  heatmap.iter(lymph_gs,
               curr.iter = i,
               samp.names = short_names,
               min.threshold = 10)
}

plot(1:number_iterations, lymph_gs[["alpha"]], type = "l",
     xlab = "Iteration", ylab = "Alpha", main = )

par(mfrow = c(4, 1))
for (i in 1:4) {
  label.switch.plot(lymph_gs, num.samples = 10)
}

dev.off()

#############################################

# Correct label switching phenomenon

lymph_ls <- correct.label.switch(lymph_gs, method.set = c("ECR"),
                                 burn.in = number_burn_in)
ls_tbl <- tibble(clust_assign = as.numeric(lymph_ls$clusters)) %>%
  mutate(mut_id = row_number())
fwrite(ls_tbl, file.path(output_dir, "clust_assign.csv"))

saveRDS(lymph_ls, file.path(output_dir, "object_lymph_ls.dat"))
save.image(file.path(output_dir, "Rsession_ls.dat"))

#############################################

# Merge clusters which differ only at level of sequencing errors
# clusters which have at least 20 muts in them, and are within 0.025 from
# each other, collapse them into a single cluster; also needs to be present
# at at least 0.025 in a sample to be considered

lymph_merge_ls <- merge.clusters.differ.by.seq.error(lymph_gs,
                                                     lymph_ls,
                                                     min_true_vaf <- 0.025,
                                                     overlap_frac <- 0.025,
                                                     min_num_muts <- 20,
                                                     ls_type <- "ECR",
                                                     number_burn_in,
                                                     version <- 2)



#############################################
ls_tbl <- tibble(clust_assign = as.numeric(lymph_merge_ls$final.cluster)) %>%
  mutate(mut_id = row_number())
ls_tbl$chrom <- lymph_depth$chrom
ls_tbl$pos <- lymph_depth$pos
fwrite(ls_tbl, file.path(output_dir, "clust_assign_posthoc.csv"))
saveRDS(lymph_merge_ls, file.path(output_dir, "object_lymph_merge_ls.dat"))

#############################################

#### for debugging at this point
if (FALSE) {
  input_dir <- "outputs/ndp_results/PD51606/"
  load(paste0(input_dir, "Rsession.dat"))
  load(paste0(input_dir, "Rsession_ls.dat"))
  lymph_merge_ls <- readRDS(paste0(input_dir, "object_lymph_merge_ls.dat"))
}

#############################################

# Final spectrum and cluster location plots
pdf(file.path(output_dir, "Cluster_and_spectrum_plots.pdf"),
    paper = "a4", width = 8, height = 10)

par_min_threshold <- 50
post_cluster_pos <- post.param.dist_posthoc(gs.out = lymph_gs,
                                            ls.merge.out = lymph_merge_ls,
                                            centiles = c(0.025, 0.5, 0.975),
                                            ls.type = "ECR",
                                            burn.in = number_burn_in,
                                            samp.names = short_names,
                                            plot.heatmap = TRUE,
                                            min.threshold = par_min_threshold)

mut_spec_by_clust <- context.extract_posthoc(mut_spec, lymph_merge_ls)

dev.off()

saveRDS(post_cluster_pos, file.path(output_dir, "object_post_cluster_pos.dat"))

## Write heatmap data to disk
cluster_tbl <- tbl_df(post_cluster_pos$heat.dat * 2)
cluster_tbl$cluster_id <- rownames(post_cluster_pos$heat.dat)

cluster_tbl$estimated.no.of.mutations <- as.numeric(
  post_cluster_pos$num.mut[
    which(post_cluster_pos$num.muts >= par_min_threshold)
  ]
)

cluster_tbl$no.of.mutations.assigned <- cluster_tbl$estimated.no.of.mutations

# Add a root node
heatmap_tbl <- rbind(cluster_tbl %>% dplyr::select(-cluster_id),
                     c(rep(1, dim(post_cluster_pos$heat.dat)[2]),
                       sum(cluster_tbl$no.of.mutations.assigned),
                       max(cluster_tbl$no.of.mutations.assigned)))

fwrite(heatmap_tbl, file.path(output_dir, "heatmap.tsv"), sep = "\t")

## Export a table with mutations per cluster
mut_per_cluster <- tibble(cluster_id = names(post_cluster_pos$num.mut),
                          num_mut = post_cluster_pos$num.mut)
fwrite(mut_per_cluster, file.path(output_dir, "muts_per_cluster.csv"))

## Export a heatmap table that includes cluster IDs
fwrite(cluster_tbl %>% dplyr::select(-estimated.no.of.mutations),
       file.path(output_dir, "cluster_and_samples.csv"))


#############################################

# Density plots of clusters with at least 100 mutations # edited to 50 muts
pdf(file.path(output_dir, "Raw_mutation_VAF_by_cluster_plots.pdf"),
    width = 12, height = 15)

which_clusters <- as.double(
  names(
    which(
      table(lymph_merge_ls[["final.cluster"]]) >= 25
    )
  )
)

for (i in which_clusters) {
  temp_plot <- cluster.density.plot_posthoc(lymph_gs,
                                            lymph_merge_ls,
                                            post_cluster_pos,
                                            i,
                                            samp.names = short_names)
  print(temp_plot)
}

dev.off()

#############################################

# Boxplot of VAFs per sample/cluster using 95.0% centiles (as conf. intervals)
# Added by Fede (fa8):
final_clusters <- unique(ls_tbl$clust_assign)
samples        <- names(post_cluster_pos$centiles[1, , 1])
num_muts       <- length(ls_tbl$clust_assign)
colores        <- rainbow(length(final_clusters))

for (sample in samples) {
  pdf(
    file.path(
      output_dir,
      sprintf("Cluster_cellular_prevalence_%s.pdf", sample)
    ),
    width = 20,
    height = 6
  )

  par(mar = c(3, 5, 1, 3))
  plot(NULL, NULL, xlim = c(0, num_muts + length(final_clusters)),
       ylim = c(0, 1.1), xaxt = "n", xlab = "", ylab = "Cell. prev")
  title(sample, line = -1)
  for (h in c(.2, .4, .6, .8)) {
    abline(h = h, col = "gray")
  }

  xpos  <- 0
  index <- 1
  for (cluster_index in final_clusters) {
    centiles <- post_cluster_pos$centiles[cluster_index, sample, ]
    centiles <- centiles * 2 # cellular prevalence
    num_muts_in_cluster <- length(which(ls_tbl$clust_assign == cluster_index))

    rect(xpos, centiles[1], xpos + num_muts_in_cluster,
         centiles[3], col = colores[index])
    lines(c(xpos, xpos + num_muts_in_cluster),
          c(centiles[2], centiles[2]), col = "black", lwd = 2)
    text((xpos + xpos + num_muts_in_cluster) / 2,
         centiles[2] + 0.075, cluster_index)
    xpos  <- xpos  + num_muts_in_cluster + 2
    index <- index + 1
  }
  dev.off()
}

save.image(file.path(output_dir, "Rsession_ls2.dat"))
