library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(tibble)
source("~/Dropbox/Single_Cell_10X/Version-2/1.1.1.cluster_mapping.R")
data_integration <- function(merge_data){
  merge_data <- lapply(X = merge_data, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  ### Select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = merge_data)
  
  ######### Perform integration #######
  immune.anchors <- FindIntegrationAnchors(object.list = merge_data, anchor.features = features)
  immune.combined <- IntegrateData(anchorset = immune.anchors)
  
  ######## Perform an integrated analysis ###########
  DefaultAssay(immune.combined) <- "integrated"
  
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindClusters(immune.combined, resolution = 0.5)
  return(immune.combined)
}
############# Visualization
data_visualization <- function(immune.combined, pdf_name){
  p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")
  p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
  p3 <- DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident",
                label = T, repel = T)
  top_row <- plot_grid(p1, p2)
  p4 <- plot_grid(top_row, p3, ncol = 1)
  pdf_real_name <- paste(pdf_name, ".pdf", sep = "")
  pdf(pdf_real_name, width = 12, height = 9) 
  plot(p4)
  dev.off()
}

############## Cluster Cell Counting #####
count_cluster_cells <- function(immune.combined, stat_name){
  merge_metadata <- immune.combined@meta.data
  metadata_stat <- merge_metadata %>% 
    group_by(orig.ident, integrated_snn_res.0.5) %>%
    dplyr::count(identified_clusters) %>%
    as.data.frame()
  
  stat_file <- paste(stat_name, "_stat.txt", sep = "")
  write.table(metadata_stat, file = stat_file, sep = "\t", row.names = F)
}

################### Seeking for marker genes ##########
get_conserved <- function(cluster){
  FindConservedMarkers(immune.combined,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       min.cells.group = 1,
                       only.pos = T) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

#####################################################################################
#####################################################################################
#####################################################################################
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/Total_Merge/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Re/HP4540Re2_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_no_doublets.RData")

merge_data <- list(HP4313Pr, HP4313Re, HP4540Pr, HP4540Re2, HP4568Pr, HP4568Re)
names(merge_data) <- c("HP4313Pr", "HP4313Re", "HP4540Pr", "HP4540Re", "HP4568Pr", "HP4568Re")
immune.combined <- data_integration(merge_data)
data_visualization(immune.combined, "Seurat_Cluster")
count_cluster_cells(immune.combined, "sepe_metadata")
DefaultAssay(immune.combined) <- "RNA"
conserved_ident_markers <- purrr::map_dfr(seq(0,22), get_conserved)
write.table(conserved_ident_markers, file = "combined_markers.txt", sep = "\t", row.names = F)

immune.combined <- RenameIdents(immune.combined,
                                "0" = "Tumor", "1" = "Tumor", "2" = "Tumor",
                                "3" = "Macrophage", "4" = "Tumor", "5" = "Tumor",
                                "6" = "Fibroblast", "7" = "T cell", "8" = "MDSC",
                                "9" = "Macrophage", "10" = "T cell", "11" = "Macrophage",
                                "12" = "Tumor", "13" = "Endothelia", "14" = "Tumor",
                                "15" = "Tumor", "16" = "B cell", "17" = "Tumor",
                                "18" = "Fibroblast", "19" = "Tumor", "20" = "T cell",
                                "21" = "Tumor", "22" = "Endothelia")
data_visualization(immune.combined, "Annot_Cluster")
immune_new_cluster_ids <- c("Tumor", "Tumor", "Tumor",
                            "Macrophage", "Tumor", "Tumor",
                            "Fibroblast", "T cell", "MDSC",
                            "Macrophage", "T cell", "Macrophage",
                            "Tumor", "Endothelia", "Tumor",
                            "Tumor", "B cell", "Tumor",
                            "Fibroblast", "Tumor", "T cell",
                            "Tumor", "Endothelia")
colnames(immune.combined@meta.data)[9] <- "identified_sep_clusters"
immune.combined <- cluster_metadata(immune.combined, immune_new_cluster_ids)
count_cluster_cells(immune.combined, "inte_metadata")
