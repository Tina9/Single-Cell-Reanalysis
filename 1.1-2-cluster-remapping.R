library(Seurat)
library(dplyr)
library(ggplot2)
options(scipen = 20)
source("~/Dropbox/Single_Cell_10X/Version-2/1.1.1.cluster_mapping.R")
############# Integrate Cluster Information #####################
label_remapping <- function(pbmc, new_cluster_ids, sample_name){
  names(new_cluster_ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new_cluster_ids)
  p1 <- DimPlot(pbmc, reduction = "umap", label = T)
  annot_name <- paste(sample_name, "_annotation.pdf", sep = "")
  pdf(annot_name, width = 9, height = 6)
  plot(p1)
  dev.off()
}

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
HP4540Pr_new_cluster_ids <- c("Tumor", "Fibroblasts", "Macrophage", "T cells",
                              "Tumor", "Tumor-proliferating", "Fibroblasts-proliferating", "Fibroblasts",
                              "Tumor", "Tumor", "Endothelia", "Undefined",
                              "B cells", "Fibroblasts", "Fibroblasts", "Neutrophils",
                              "Macrophage")
label_remapping(HP4540Pr, HP4540Pr_new_cluster_ids, "HP4540Pr")
HP4540Pr <- cluster_metadata(HP4540Pr, HP4540Pr_new_cluster_ids)
save(HP4540Pr, file = "HP4540Pr_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Re/HP4540Re2_no_doublets.RData")
HP4540Re_new_cluster_ids <- c("Tumor", "Tumor", "Tumor", "Macrophage",
                              "Tumor", "T cells", "Tumor", "Tumor",
                              "Tumor-proliferating", "Endothelia", "Fibroblasts",
                              "Tumor", "stromal", "T cells/NKT", "Epithelia",
                              "Neutrophils", "stromal")
label_remapping(HP4540Re2, HP4540Re_new_cluster_ids, "HP4540Re2")
HP4540Re2 <- cluster_metadata(HP4540Re2, HP4540Re_new_cluster_ids)
save(HP4540Re2, file = "HP4540Re2_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_no_doublets.RData")
HP4313Pr_new_cluster_ids <- c("Tumor", "Tumor-proliferating", "Tumor", "Macrophage",
                              "Tumor", "Tumor", "MDSCs", "T cells",
                              "Macrophage", "T cells", "Fibroblasts", "Fibroblasts",
                              "stromal", "Endothelia", "DC")
label_remapping(HP4313Pr, HP4313Pr_new_cluster_ids, "HP4313Pr")
HP4313Pr <- cluster_metadata(HP4313Pr, HP4313Pr_new_cluster_ids)
save(HP4313Pr, file = "HP4313Pr_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_no_doublets.RData")
HP4313Re_new_cluster_ids <- c("Tumor", "Tumor", "Macrophage", "Tumor",
                              "Tumor-proliferating", "Tumor", "MDSCs",
                              "T cells", "Macrophage", "Fibroblasts", "Undefined",
                              "Endothelia")
label_remapping(HP4313Re, HP4313Re_new_cluster_ids, "HP4313Re")
HP4313Re <- cluster_metadata(HP4313Re, HP4313Re_new_cluster_ids)
save(HP4313Re, file = "HP4313Re_no_doublets.RData")


setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_no_doublets.RData")
HP4568Pr_new_cluster_ids <- c("Undefined", "Tumor", "Tumor", "Tumor",
                              "Tumor", "Endothelia", "Macrophage",
                              "Tumor", "Undefined", "Tumor", "Undefined",
                              "Tumor", "Fibroblasts", "T cells", "MDSC/Neutrophils",
                              "Endothelia")
label_remapping(HP4568Pr, HP4568Pr_new_cluster_ids, "HP4568Pr")
HP4568Pr <- cluster_metadata(HP4568Pr, HP4568Pr_new_cluster_ids)
save(HP4568Pr, file = "HP4568Pr_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_no_doublets.RData")
HP4568Re_new_cluster_ids <- c("Tumor", "Macrophage", "Tumor", "Tumor", "Tumor-proliferating",
                              "T cells", "Neutrophils", "Tumor", "B cells", "T cells",
                              "Endothelia", "Macrophage", "Fibroblasts", "Macrophage-proliferating",
                              "Fibroblasts")
label_remapping(HP4568Re, HP4568Re_new_cluster_ids, "HP4568Re")
HP4568Re <- cluster_metadata(HP4568Re, HP4568Re_new_cluster_ids)
save(HP4568Re, file = "HP4568Re_no_doublets.RData")

