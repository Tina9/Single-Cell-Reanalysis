####### Integrating Cluster Information into Metadata #################
cluster_metadata <- function(pbmc, identified_clusters){
  seurat_clusters <- levels(pbmc@meta.data$seurat_clusters)
  cluster_refer <- as.data.frame(cbind(seurat_clusters, identified_clusters))
  meta_data <- pbmc@meta.data
  meta_data$identified_clusters <- cluster_refer[match(meta_data$seurat_clusters, cluster_refer$seurat_clusters), "identified_clusters"]
  pbmc@meta.data <- meta_data
  return(pbmc)
}


# load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
# pbmc <- HP4540Pr
# identified_clusters <- c("Tumor", "Fibroblasts", "Macrophage", "T cells",
#                          "Tumor", "Tumor", "Fibroblasts", "Fibroblasts",
#                          "Tumor", "Tumor", "Endothelia", "Undefined",
#                          "B cells", "Fibroblasts", "Fibroblasts", "MDSCs",
#                          "Macrophage")
# pbmc <- cluster_metadata(pbmc, identified_clusters)

