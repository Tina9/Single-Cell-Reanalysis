library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(SingleR)
options(scipen = 20)

read_data <- function(data_dir, sample_name, doublets_path){
  pbmc.data <- Read10X(data_dir)
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = sample_name)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt")
  plot_data <- pbmc@meta.data
  
  # Using density plot to decide the conditions
  p1 <- ggplot(plot_data, aes(x=nCount_RNA, color=orig.ident)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none")
  p2 <- ggplot(plot_data, aes(x = nFeature_RNA, color = orig.ident)) + 
    geom_density() +
    theme_bw() +
    theme(legend.position = "none")
  p3 <- ggplot(plot_data, aes(x = percent.mt, color = orig.ident)) +
    geom_density() +
    theme_bw()
  
  filter_plot_name <- paste(sample_name, "_filter_density_plot.pdf", sep = "")
  pdf(filter_plot_name, height = 4, width = 12)
  plot(p1 + p2 + p3)
  dev.off()
  
  ## Integrate the information of doublets into the metadata
  doublets <- read.table(doublets_path, header = T, sep = ",", row.names = 1)
  pbmc@meta.data$doublets <- doublets$predicted_doublets
  pbmc@meta.data$doublets_score <- doublets$doublet_scores
  
  return(pbmc)
}

# Using doublets as the condition
doublets_removing <- function(pbmc, min_feature, max_feature, mt_percent){
  pbmc <- subset(pbmc, subset = nFeature_RNA > min_feature &
                   nFeature_RNA < max_feature &
                   percent.mt < mt_percent &
                   doublets == "False")
  return(pbmc)
}

# Do not use doublets as the condition
doublets_no_removing <- function(pbmc, min_feature, max_feature, mt_percent){
  pbmc <- subset(pbmc, subset = nFeature_RNA > min_feature & 
                   nFeature_RNA < max_feature &
                   percent.mt < mt_percent)
  return(pbmc)
}

data_correction <- function(pbmc, sample_name){
  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, dims = 1:35)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  pbmc <- RunUMAP(pbmc, dims = 1:35)
  p1 <- DimPlot(pbmc, reduction = "umap", label = T)
  p2 <- FeaturePlot(pbmc, features = "nFeature_RNA", cols = c("white", "red"))
  cluster_umap <- paste(sample_name, "_cluster_umap.pdf", sep = "")
  pdf(cluster_umap, height = 6, width = 16)
  plot(p1 + p2)
  dev.off()
  return(pbmc)
}

markers_seeking <- function(pbmc, sample_name){
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
  markers <- pbmc.markers %>% group_by(cluster) %>% as.data.frame
  mark_file <- paste(sample_name, "_markers.txt", sep = "")
  write.table(markers, file = mark_file, sep = "\t", row.names = F)
  
  return(pbmc)
}

#####################################################################
#####################################################################
#####################################################################
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/")
dir_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4540Pr/"
doublets_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4540Pr/doublet.txt"
## Removing doublets
HP4540Pr <- read_data(dir_path, "HP4540Pr", doublets_path)
HP4540Pr <- doublets_removing(HP4540Pr, 200, 7000, 25)
HP4540Pr <- data_correction(HP4540Pr, "HP4540Pr_no_doublets")
HP4540Pr <- markers_seeking(HP4540Pr, "HP4540Pr_no_doublets")
save(HP4540Pr, file = "HP4540Pr_no_doublets.RData")
## Do not remove doublets
HP4540Pr <- read_data(dir_path, "HP4540Pr", doublets_path)
HP4540Pr <- doublets_no_removing(HP4540Pr, 200, 7000, 25)
HP4540Pr <- data_correction(HP4540Pr, "HP4540Pr_with_doublets")
HP4540Pr <- markers_seeking(HP4540Pr, "HP4540Pr_with_doublets")
save(HP4540Pr, file = "HP4540Pr_with_doublets.RData")

#####################################################################
#####################################################################
#####################################################################
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/")
dir_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Pr/"
doublets_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Pr/doublet.txt"
## Removing doublets
HP4313Pr <- read_data(dir_path, "HP4313Pr", doublets_path)
HP4313Pr <- doublets_removing(HP4313Pr, 200, 9500, 25)
HP4313Pr <- data_correction(HP4313Pr, "HP4313Pr_no_doublets")
HP4313Pr <- markers_seeking(HP4313Pr, "HP4313Pr_no_doublets")
save(HP4313Pr, file = "HP4313Pr_no_doublets.RData")
## Do not remove doublets
HP4313Pr <- read_data(dir_path, "HP4313Pr", doublets_path)
HP4313Pr <- doublets_no_removing(HP4313Pr, 200, 9500, 25)
HP4313Pr <- data_correction(HP4313Pr, "HP4313Pr_with_doublets")
HP4313Pr <- markers_seeking(HP4313Pr, "HP4313Pr_with_doublets")
save(HP4313Pr, file = "HP4313Pr_with_doublets.RData")

#####################################################################
#####################################################################
#####################################################################
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/")
dir_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Re/"
doublets_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4313Re/doublet.txt"
## Removing doublets
HP4313Re <- read_data(dir_path, "HP4313Re", doublets_path)
HP4313Re <- doublets_removing(HP4313Re, 200, 9000, 25)
HP4313Re <- data_correction(HP4313Re, "HP4313Re_no_doublets")
HP4313Re <- markers_seeking(HP4313Re, "HP4313Re_no_doublets")
save(HP4313Re, file = "HP4313Re_no_doublets.RData")
## Do not remove doublets
HP4313Re <- read_data(dir_path, "HP4313Re", doublets_path)
HP4313Re <- doublets_no_removing(HP4313Re, 200, 9000, 25)
HP4313Re <- data_correction(HP4313Re, "HP4313Re_with_doublets")
HP4313Re <- markers_seeking(HP4313Re, "HP4313Re_with_doublets")
save(HP4313Re, file = "HP4313Re_with_doublets.RData")

#####################################################################
#####################################################################
#####################################################################
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/")
dir_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Pr/"
doublets_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Pr/doublet.txt"
## Removing doublets
HP4568Pr <- read_data(dir_path, "HP4568Pr", doublets_path)
HP4568Pr <- doublets_removing(HP4568Pr, 200, 9500, 25)
HP4568Pr <- data_correction(HP4568Pr, "HP4568Pr_no_doublets")
HP4568Pr <- markers_seeking(HP4568Pr, "HP4568Pr_no_doublets")
save(HP4568Pr, file = "HP4568Pr_no_doublets.RData")
## Do not remove doublets
HP4568Pr <- read_data(dir_path, "HP4568Pr", doublets_path)
HP4568Pr <- doublets_no_removing(HP4568Pr, 200, 9500, 25)
HP4568Pr <- data_correction(HP4568Pr, "HP4568Pr_with_doublets")
HP4568Pr <- markers_seeking(HP4568Pr, "HP4568Pr_with_doublets")
save(HP4568Pr, file = "HP4568Pr_with_doublets.RData")
#####################################################################
#####################################################################
#####################################################################
setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/")
dir_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Re/"
doublets_path <- "~/Dropbox/BRCA1-PARPi-10X/Seurat_Raw_Data/HP4568Re/doublet.txt"
## Removing doublets
HP4568Re <- read_data(dir_path, "HP4568Re", doublets_path)
HP4568Re <- doublets_removing(HP4568Re, 200, 8500, 25)
HP4568Re <- data_correction(HP4568Re, "HP4568Re_no_doublets")
HP4568Re <- markers_seeking(HP4568Re, "HP4568Re_no_doublets")
save(HP4568Re, file = "HP4568Re_no_doublets.RData")
## Do not remove doublets
HP4568Re <- read_data(dir_path, "HP4568Re", doublets_path)
HP4568Re <- doublets_no_removing(HP4568Re, 200, 8500, 25)
HP4568Re <- data_correction(HP4568Re, "HP4568Re_with_doublets")
HP4568Re <- markers_seeking(HP4568Re, "HP4568Re_with_doublets")
save(HP4568Re, file = "HP4568Re_with_doublets.RData")
