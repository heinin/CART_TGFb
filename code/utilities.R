lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

nullToNA <- function(x) {
  x[sapply(x, is.null)] <- NA
  return(x)
}

get_pcs <- function(seurat_object, reduction_name="pca") {
  
  # Determine percent of variation associated with each PC
  pct <- seurat_object[[reduction_name]]@stdev / sum(seurat_object[[reduction_name]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % 
  # variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co1
  
  # Determine the difference between variation of PC and subsequent PC and
  # selecting last point where change of % of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Minimum of the two calculation
  #pcs <- min(co1, co2)
  
  c(co1, co2)
  
}

# TODO: move to the function, save to a file
# Create a dataframe with values
#plot_df <- data.frame(pct = pct, 
#                      cumu = cumu, 
#                      rank = 1:length(pct))
#
## Elbow plot to visualize 
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()

potato_pcs <- function(seurat_object, filename, pcs){
  # Empty matrices for the distances for the UMAP cell distances
  potato <-  matrix(ncol = length(pcs), nrow = ncol(seurat_object))
  potato2 <-  matrix(ncol = length(pcs), nrow = ncol(seurat_object))

  j <- 0
  for(i in pcs) {
    j <- j + 1
    print(i)
    #print(j)
    # Running UMAP with a given set of PCs
    onion <- RunUMAP(seurat_object, dims = 1:i, verbose = T)
    # UMAP distances between cells are stored in cell.embeddings
    potato[, j] <- onion@reductions$umap@cell.embeddings[, 1]
    potato2[, j] <- onion@reductions$umap@cell.embeddings[, 2]
  }
  
  colnames(potato) <- pcs
  colnames(potato2) <- pcs
  
  # Subsetting to a random set of cells
  mrand_obs <- NULL
  sample_subset <- sample(1:ncol(seurat_object), 3000)
  #potato_new <- potato[sample(1:ncol(seurat_object), 3000), ]
  #potato2_new <- potato2[sample(1:ncol(seurat_object), 3000), ]
  potato_new <- potato[sample_subset, ]
  potato2_new <- potato2[sample_subset, ]
  
  # Mantel's test (correlation between two distance matrices)
  # Comparing the distances between neighboring UMAPs
  for(i in  1:c(dim(potato_new)[2] - 1)) {
    onion <- dist(potato_new[, i])
    onion2 <- dist(potato_new[, i + 1])
    onion3 <- mantel.randtest(onion, onion2)$obs
    onion <- dist(potato2_new[, i])
    onion2 <- dist(potato2_new[, i + 1])
    onion4 <- mantel.randtest(onion, onion2)$obs
    
    #plot(r1 <- mantel.randtest(onion,onion2), main = "Mantel's test")
    #r1
    
    mrand_obs <- rbind(mrand_obs, cbind(onion3, onion4))
  }
  
  pdf(filename)
  #ylim = c(0.8, 1)
  plot(pcs[1:length(pcs)-1],
       mrand_obs[, 1], xlab="PCs included in the UMAP", ylab="Mantel's test observed correlation",
       ylim = c(0.2, 1.2))
  lines(spline(x = pcs[1:length(pcs)-1],
               y = mrand_obs[, 1]))
  # ylim = c(0.8, 1)
  points(pcs[1:length(pcs)-1],
         mrand_obs[, 2], col = "red", ylim = c(0.2, 1.2))
  lines(spline(x = pcs[1:length(pcs)-1],
               y = mrand_obs[, 2]), col = "red")
  dev.off()
  
}