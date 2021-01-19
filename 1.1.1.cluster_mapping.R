####### Integrating Cluster Information into Metadata #################
cluster_metadata <- function(pbmc, identified_clusters){
  seurat_clusters <- levels(pbmc@meta.data$seurat_clusters)
  cluster_refer <- cbind(seurat_clusters, identified_clusters)
  cluster_meta <- merge(pbmc@meta.data, cluster_refer, 
                        by.x = "seurat_clusters",
                        by.y = "seurat_clusters",
                        all.x = T)
  pbmc@meta.data <- cluster_meta
  return(pbmc)
}