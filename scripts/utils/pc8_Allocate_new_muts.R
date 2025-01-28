  # author: Peter Campbell 
  # Function to allocate a new set of mutations to an existing DP clustering based on VAFs
  # Only allocates to non-zero clusters
  # Uses a naive Bayesian classifier - weighted by number of muts in cluster and VAF means
  # In this sense, represents a true posterior classification
  
  # ls.merge.out is the output from the DP clustering after merge.clusters.differ.by.seq.error() function
  # post.ls.dist is the posterior distribution of cluster positions from post.param.dist()
  # new.y and new.N are matrices of new mutations to be allocated to the best clusters. 
  # The columns are in the same order as the original y and N for the DP (ie samples are in same order)
  # new.y has number of reads reporting the variant; new.N has total depth
  
  # Returns vector of best allocations for the new mutations

allocate.new.muts.to.cluster <- function(ls.merge.out, post.ls.dist, new.y, new.N) {  
  # Extract variables and data
  ls.orig <- ls.merge.out[["final.cluster"]]
  post.pi <- post.ls.dist[["centiles"]][,,"Centile.0.5"]
  post.pi <- t(post.pi)
  new.N.minus.y <- new.N - new.y
  num.new.muts <- nrow(new.y)
  
  # Generate weights for each cluster
  cluster.counts <- table(ls.orig)
  log.weights <- log(cluster.counts)
  num.clusters <- length(cluster.counts)
  
  # Extract mean VAFs for each non-zero cluster
  post.pi <- post.pi[, colnames(post.pi) %in% paste0("Cl.",names(cluster.counts))]
  log.post.pi <- log(post.pi)
  log.1.minus.pi <- log(1 - post.pi)
  
  # Allocate new mutations to clusters
  y.ij.log.post.pi <- new.y %*% log.post.pi
  N.minus.y.log.1.minus.pi <- new.N.minus.y %*% log.1.minus.pi
  log.Pr.S <- matrix(rep(log.weights, each=num.new.muts), nrow=num.new.muts, byrow = FALSE) +
    y.ij.log.post.pi +
    N.minus.y.log.1.minus.pi
  Pr.S <- exp(log.Pr.S) 
  Pr.S[rowSums(Pr.S) == 0,] <- 1/num.clusters
  
  S.i <- sapply(1:num.new.muts, function(k) {sum(rmultinom(1,1,Pr.S[k,]) * (1:num.clusters))})
  S.i <- as.double(names(cluster.counts)[S.i])
  return(S.i)
}

