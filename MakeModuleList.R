# mat = matrix of normalized RNA-seq expression data where each row is a gene and each column is a cluster/cell-type.
# Cluster data derived from single cells can be median gene expression across all cells within the cluster or sum counts subsampled for a fixed number of cells (e.g. 1000 cells).
# sd_cutoff = number of standard deviations above the mean expression across all clusters/celltypes for which you want to specify a gene as being included in a module.
# min_genes = minimum genes that are needed to be in a module for a module to be defined

MakeModList <- function(mat, sd_cutoff = 2, min_genes = 10){
  ScaledMatrix <- scale(t(mat))
  foo <- c()
  for(i in rownames(ScaledMatrix)){
    foo[[i]] <- ScaledMatrix[i,]
    foo[[i]] <- foo[[i]][foo[[i]] > sd_cutoff  & !is.na(foo[[i]])]
    foo[[i]] <- foo[[i]][order(foo[[i]], decreasing = TRUE)]
    if(length(foo[[i]]) < min_genes){
      foo[[i]] <- NULL
      }
    foo[[i]] <- names(foo[[i]])
  }
  module_list <- foo
  print(module_list)
}

