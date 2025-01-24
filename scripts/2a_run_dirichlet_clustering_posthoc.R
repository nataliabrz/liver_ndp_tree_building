# Run Dirichlet clustering
# --
# Runs Dirichlet clustering and label switching on low input mutation data
#
# Args (command line):
#   1: Path to git repository (ie. parent directory where repository was cloned into)
#   2: Path to CSV containing counts of alt base (rows are variants, any columns starting with 'PD' are samples)
#   3: Path to CSV containing depth (rows are variants, any columns starting with 'PD' are samples)
#   4: Path to table with mutation contexts (generated with Peter Campbell's "Context_pull_build37.pl" script)
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
args = commandArgs(trailingOnly=TRUE)
args

# Check for command line args, else set defaults
if (length(args)==0) {
  repo_location = "scripts/utils/"
  fpath_alt_csv = "data/filtered_calls/snv/PD51606_ndp_alt_bb_flt.csv"
  fpath_depth_csv = "data/filtered_calls/snv/PD51606_ndp_depth_bb_flt.csv"
  fpath_mutant_locations = "data/filtered_calls/snv/PD51606_mut_context_GRCh38.txt"
  sex.file = "data/metadata/sex_file.txt"
  output_dir = "outputs/ndp_results/PD51606"
  num_burnin = 15000
} else if (length(args)==7) {
  repo_location = args[1]
  fpath_alt_csv = args[2]
  fpath_depth_csv = args[3]
  fpath_mutant_locations = args[4]
  sex.file = args[5]
  output_dir = args[6]
  num_burnin = as.numeric(args[7])
} else {
  stop("Script takes either exactly 0 or 7 arguments.", call.=FALSE)
}

# Report arguments
message(sprintf('Reading alt counts from %s', fpath_alt_csv))
message(sprintf('Reading depth from %s', fpath_depth_csv))
message(sprintf('Reading mutation contexts from %s', fpath_mutant_locations))
message(sprintf('Writing output to %s', output_dir))
message(sprintf('Number of burnin cycles is %s', num_burnin))

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


source(file.path(repo_location,'Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R'))
source(file.path(repo_location,'Mutation_Cluster_Dirichlet_posthoc_merge_fx.R'))

dir.create(output_dir, showWarnings = F)

# Define params
number.iterations <- 25000
number.burn.in <- 15000
number.burn.in <- num_burnin

set.seed(28)

## define gender
f_list <- c("F", "f", "FEMALE", "female", "Female")
m_list<- c("M", "m", "MALE", "male", "Male")

#############################################
## -- Run clustering --

lymph.var <- read.table(fpath_alt_csv, sep=",", header=TRUE, stringsAsFactors = FALSE)

lymph.depth <- read.table(fpath_depth_csv, sep=",", header=TRUE, stringsAsFactors = FALSE)

mut.spec <- read.table(fpath_mutant_locations, sep="\t", header=TRUE, 
                       stringsAsFactors = FALSE)

######### added by nb15:
#for males, double the coverage on X and Y
#for females, remove chr Y muts
sex_tbl <- fread(sex.file)
this_patient=substr(colnames(lymph.depth)[grep("PD",colnames(lymph.depth))][1],1,7)
pat_sex = sex_tbl[sex_tbl$patientID == this_patient]$sex
if (pat_sex %in% m_list) {
cat("patient is male, doubling coverage on X and Y chromosomes")
lymph.depth[lymph.depth$chrom %in% c("chrX","chrY"),grep("PD",colnames(lymph.depth))] <- lymph.depth[lymph.depth$chrom %in% c("chrX","chrY"),grep("PD",colnames(lymph.depth))]*2
}else if (pat_sex %in% f_list){
cat("patient is female, removing muts on Y chromosome")
lymph.var <- lymph.var[lymph.var$chrom != "chrY",] # Remove Y chrom muts
lymph.depth <- lymph.depth[lymph.depth$chrom != "chrY",]
mut.spec <- mut.spec[ mut.spec$chrom != "chrY",]
}

#########

lymph.y <- as.matrix(lymph.var[,which(grepl('PD', names(lymph.var)))])
lymph.N <- as.matrix(lymph.depth[,which(grepl('PD', names(lymph.depth)))])
lymph.gs <- muts.cluster.dirichlet.gibbs(C = 100, y = lymph.y, N = lymph.N, iter=number.iterations,
                                         iter.chunk=100, outdir=output_dir)

short.names <- names(lymph.var)[which(grepl('PD', names(lymph.var)))]
short.names[2:length(short.names)] <- substring(short.names[2:length(short.names)], 8)

save.image(file.path(output_dir, 'Rsession.dat'))

#############################################

# Generate convergence plots

pdf(file.path(output_dir, "Convergence_plots.pdf"), paper="a4", width = 8, height = 10)

for (i in c(10, 50, 100, 200, 500, seq(1000, number.iterations, 1000))) {
  heatmap.iter(lymph.gs, curr.iter = i, samp.names = short.names, min.threshold = 10)
}

plot(1:number.iterations, lymph.gs[["alpha"]], type="l", xlab="Iteration", ylab = "Alpha", main = )

par(mfrow = c(4,1))
for (i in 1:4) {
  label.switch.plot(lymph.gs, num.samples = 10)
}

dev.off()

#############################################

# Correct label switching phenomenon

lymph.ls <- correct.label.switch(lymph.gs, method.set = c("ECR"), burn.in = number.burn.in)
ls_tbl = tibble(clust_assign = as.numeric(lymph.ls$clusters)) %>%
  mutate(mut_id = row_number())
fwrite(ls_tbl, file.path(output_dir, 'clust_assign.csv'))

saveRDS(lymph.ls, file.path(output_dir, 'object_lymph.ls.dat'))
save.image(file.path(output_dir, 'Rsession_ls.dat'))

#############################################

# Merge clusters which differ only at level of sequencing errors
# clusters which have at least 20 muts in them, and are within 0.025 from each other, collapse them into a single cluster; also needs to be present at at least 0.025 in a sample to be considered
#lymph.merge.ls <- merge.clusters.differ.by.seq.error(lymph.gs, lymph.ls, min.true.VAF = 0.025, overlap.frac=0.025, min.num.muts = 20, ls.type="ECR", number.burn.in, version=2)
## nb15 I am printing the function here because I am getting different results when running as function vs. line by line... this is a temporary fix...
gs.out=lymph.gs
ls.out=lymph.ls
min.true.VAF = 0.025
overlap.frac=0.025
min.num.muts = 20
ls.type="ECR"
burn.in=number.burn.in
version=2

ls.final <- ls.out$clusters[1,]
  unique.clusters <- table(ls.final)
  kept_clusters = which(unique.clusters >= min.num.muts)
  unique.clusters <- unique.clusters[kept_clusters]
  unique.clusters <- names(unique.clusters)[order(unique.clusters, decreasing = TRUE)]

  y <- gs.out[["y1"]]
  N <- gs.out[["N1"]]
  merge.spl <- gs.out[["merge.split"]]

  ls.perm <- ls.out$permutations[ls.type][[1]]
  pi.h.j <- gs.out[["pi.h.j"]]
  num.clusters <- dim(pi.h.j)[3]
  num.samples <- dim(pi.h.j)[2]
  num.iter <- dim(pi.h.j)[1]

  pi.mcmc <- array(NA, dim = c(num.iter - burn.in, num.clusters, num.samples))
  for (i in (burn.in+1):num.iter) {pi.mcmc[i-burn.in,,] <- t(pi.h.j[i,,])}
  # Permute according to label-switching() defined optimal labelling 
  pi.mcmc <- permute.mcmc(mcmc = pi.mcmc, permutations = ls.perm)$output

  # For successful merge-split steps after burn.in, turn estimates of pi before the merge-split to NA
  merge.spl <- merge.spl[(burn.in+1):num.iter,]
  merge.spl <- merge.spl[merge.spl$decision,]

  # condition to handle cases where no merging of clusters is required
  if(dim(merge.spl)[1] > 0){ #nb15 commented out, I have cases where clusters should be merged regardless

  merge.spl$new.cluster.loc <- as.numeric(merge.spl$new.cluster.loc)
  for (i in 1:nrow(merge.spl)) {
    final.cluster.1 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$first.cluster[i]]
    final.cluster.2 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$second.cluster[i]]
    final.cluster.3 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$new.cluster.loc[i]]

    # instead of assigning N/A's, assign 0's instead to ensure no N/A's are present
    # during subsequent mean.vaf calculation and heatmap plotting ???
    if (merge.spl$first.cluster[i] == merge.spl$second.cluster[i]) {
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.1, ] <- NA
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.3, ] <- NA
    } else {
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.1, ] <- NA
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.2, ] <- NA
    }
  }

  }

  # since above can assign some N/A's to pi.mcmc some columns in mean.vaf will be N/A,
  # should N/A's be converted to zeros, or should merge step be skipped if N/A's detected ???
  mean.vaf <- apply(pi.mcmc, MARGIN = 2, FUN = colMeans, na.rm=TRUE)
  if (version == 2) {
    S.i <- read.table(gs.out[["S.i.filename"]], header=TRUE, sep="\t", stringsAsFactors = FALSE,
                      skip = burn.in)
    size.clusters <- t(apply(S.i, MARGIN = 1,
                             FUN = function(a) {table(factor(a, levels = 1:num.clusters))}))
  }

for (i in 1:(length(unique.clusters)-1)) {
    ind.i <- as.double(unique.clusters[i])
    # to handle the case of N/A values in mean.vaf columns:
    # skip merge step if N/A values detected
    if(sum(is.na(mean.vaf[,ind.i])) > 0){next}

    for (j in (i+1):length(unique.clusters)) {
      ind.j <- as.double(unique.clusters[j])
      # to handle the case of N/A values in mean.vaf columns:
      # skip merge step if N/A values detected
      if(sum(is.na(mean.vaf[,ind.j])) > 0){next}
      # zero out NA's
      # if(length(is.na(mean.vaf[,ind.i])) > 0){
      #  idx = which(is.na(mean.vaf[,ind.i]))
      #  mean.vaf[idx,ind.i] = 0
      # }

      nonzero.samps <- which(mean.vaf[,ind.i] >= min.true.VAF | mean.vaf[,ind.j] >= min.true.VAF )
      num.overlaps <- sapply(1:dim(pi.mcmc)[3], function(a,k,x,y) {
        temp <- na.exclude(a[,x,k] > a[,y,k]);
        if (length(temp)>0) {sum(temp) / length(temp)} else {0}},
        a=pi.mcmc, x=ind.i, y=ind.j)
      if (version == 1 &
          all(num.overlaps[nonzero.samps] > overlap.frac &
              num.overlaps[nonzero.samps] < 1-overlap.frac)) {
        # Clusters should be merged
        # Merge them into the larger cluster (ind.i)
        ls.final[ls.final == ind.j] <- ind.i
      } else if (version == 2 & all(abs(mean.vaf[,ind.i] - mean.vaf[,ind.j]) < overlap.frac)) {
        print(paste0("cluster 1 is ",ind.i))
        print(paste0("cluster 2 is ",ind.j))
        print(paste0("clusters not differing by ",overlap.frac," in any sample, merging clusters"))
        # Clusters should be merged - weight the pi estimates by size of cluster
        ls.final[ls.final == ind.j] <- ind.i
        pi.mcmc[ , ind.i, ] <- (size.clusters[, ind.i] * pi.mcmc[ , ind.i, ] +
                                  size.clusters[, ind.j] * pi.mcmc[ , ind.j, ]) /
          (size.clusters[, ind.i] + size.clusters[, ind.j])
      print(paste0(length(unique(ls.final))," clusters after updating"))
        }
    }
  }

lymph.merge.ls <- list(final.cluster = ls.final, pi.mcmc = pi.mcmc)
#############################################
ls_tbl = tibble(clust_assign = as.numeric(lymph.merge.ls$final.cluster)) %>%
  mutate(mut_id = row_number())
ls_tbl$chrom = lymph.depth$chrom
ls_tbl$pos = lymph.depth$pos
fwrite(ls_tbl, file.path(output_dir, 'clust_assign_posthoc.csv'))
saveRDS(lymph.merge.ls, file.path(output_dir, 'object_lymph.merge.ls.dat'))

#############################################

#### for debugging at this point
if(FALSE){
input_dir="outputs/ndp_results/PD51606/"
load(paste0(input_dir,"Rsession.dat"))
load(paste0(input_dir,"Rsession_ls.dat"))
lymph.merge.ls <- readRDS(paste0(input_dir,"object_lymph.merge.ls.dat"))
}

#############################################

# Final spectrum and cluster location plots
pdf(file.path(output_dir, "Cluster_and_spectrum_plots.pdf"), paper="a4", width = 8, height=10)

par_min.threshold = 50
post.cluster.pos <- post.param.dist_posthoc(gs.out = lymph.gs, ls.merge.out = lymph.merge.ls, 
                                    centiles = c(0.025, 0.5, 0.975), ls.type = "ECR", 
                                    burn.in = number.burn.in, samp.names = short.names, 
                                    plot.heatmap = TRUE, min.threshold = par_min.threshold)  

mut.spec.by.clust <- context.extract_posthoc(mut.spec, lymph.merge.ls)

dev.off()

saveRDS(post.cluster.pos, file.path(output_dir, 'object_post.cluster.pos.dat'))

## Write heatmap data to disk
cluster_tbl = tbl_df(post.cluster.pos$heat.dat*2)
cluster_tbl$cluster_id = rownames(post.cluster.pos$heat.dat)

# to handle case of NA clusters in post.param.dist_posthoc
# but NA clusters should disappear by increasing iterations to 16k
# est.num.muts = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])
# if(length(post.cluster.pos$idx.na) > 0){
#  est.num.muts = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])[-post.cluster.pos$idx.na]
# }

cluster_tbl$estimated.no.of.mutations = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])
# cluster_tbl$estimated.no.of.mutations = est.num.muts
cluster_tbl$no.of.mutations.assigned = cluster_tbl$estimated.no.of.mutations

# Add a root node 
heatmap_tbl = rbind(cluster_tbl %>% dplyr::select(-cluster_id),
                    c(rep(1, dim(post.cluster.pos$heat.dat)[2]), sum(cluster_tbl$no.of.mutations.assigned), max(cluster_tbl$no.of.mutations.assigned)))

fwrite(heatmap_tbl, file.path(output_dir, 'heatmap.tsv'), sep='\t')

## Export a table with mutations per cluster
mut_per_cluster = tibble(cluster_id = names(post.cluster.pos$num.mut), num_mut = post.cluster.pos$num.mut)
fwrite(mut_per_cluster, file.path(output_dir, 'muts_per_cluster.csv'))

## Export a heatmap table that includes cluster IDs
fwrite(cluster_tbl %>% dplyr::select(-estimated.no.of.mutations), file.path(output_dir, 'cluster_and_samples.csv'))


#############################################

# Density plots of clusters with at least 100 mutations # edited to 50 muts
pdf(file.path(output_dir, "Raw_mutation_VAF_by_cluster_plots.pdf"), width = 12, height = 15)

which.clusters <- as.double(names(which(table(lymph.merge.ls[["final.cluster"]]) >= 25)))

for (i in which.clusters) {
  temp.plot <- cluster.density.plot_posthoc(lymph.gs, lymph.merge.ls, post.cluster.pos, i, samp.names = short.names)
  print(temp.plot)
}

dev.off()

#############################################

# Boxplot of VAFs per sample/cluster using 95.0% centiles (as conf. intervals)
# Added by Fede (fa8):
final_clusters <- unique(ls_tbl$clust_assign)
samples        <- names(post.cluster.pos$centiles[1,,1])
num_muts       <- length(ls_tbl$clust_assign)
colores        <- rainbow(length(final_clusters))

for(sample in samples) {
  pdf(file.path(output_dir, sprintf("Cluster_cellular_prevalence_%s.pdf", sample)), width = 20, height=6)
  #par(mfrow=c(length(samples),1))
  par(mar=c(3,5,1,3))
  plot(NULL,NULL,xlim=c(0,num_muts+length(final_clusters)),ylim=c(0,1.1),xaxt='n',xlab="",ylab="Cell. prev")
  title(sample, line = -1)
  for(h in c(.2,.4,.6,.8)) {
    abline(h=h,col="gray");
  }
  xpos  <- 0
  index <- 1
  for(cluster_index in final_clusters) {
    centiles <- post.cluster.pos$centiles[cluster_index,sample,]
    centiles <- centiles * 2 # cellular prevalence
    num_muts_in_cluster <- length(which(ls_tbl$clust_assign==cluster_index))
    #rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
    #     col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
    #     ...)
    rect(xpos,centiles[1],xpos+num_muts_in_cluster,centiles[3],col=colores[index])
    lines(c(xpos,xpos+num_muts_in_cluster),c(centiles[2],centiles[2]),col="black",lwd=2)
    text((xpos + xpos + num_muts_in_cluster)/2,centiles[2]+0.075,cluster_index)
    xpos  <- xpos  + num_muts_in_cluster + 2
    index <- index + 1
  }
  dev.off()
}


#save.image()
save.image(file.path(output_dir, 'Rsession_ls2.dat'))
