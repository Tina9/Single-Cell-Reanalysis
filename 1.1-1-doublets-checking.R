library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
options(scipen = 20)

doublets_checking <- function(pbmc_no_doublets, pbmc_with_doublets, sample_name){
  p1 <- DimPlot(pbmc_no_doublets, reduction = "umap", label = T)
  p2 <- FeaturePlot(pbmc_no_doublets, features = "doublets_score", cols = c("white", "red"))
  p3 <- DimPlot(pbmc_with_doublets, reduction = "umap", label = T)
  p4 <- FeaturePlot(pbmc_with_doublets, features = "doublets_score", cols = c("white", "red"))
  p5 <- plot_grid(p1, p2, p3, p4)
  
  doublets_pdf <- paste(sample_name, "_doublets.pdf", sep = "")
  pdf(doublets_pdf, width = 12, height = 9)
  plot(p5)
  dev.off()
}

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
HP4540Pr_no_doublets <- HP4540Pr
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_with_doublets.RData")
HP4540Pr_with_doublets <- HP4540Pr
doublets_checking(HP4540Pr_no_doublets, HP4540Pr_with_doublets, "HP4540Pr")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_no_doublets.RData")
HP4313Pr_no_doublets <- HP4313Pr
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_with_doublets.RData")
HP4313Pr_with_doublets <- HP4313Pr
doublets_checking(HP4313Pr_no_doublets, HP4313Pr_with_doublets, "HP4313Pr")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_no_doublets.RData")
HP4313Re_no_doublets <- HP4313Re
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_with_doublets.RData")
HP4313Re_with_doublets <- HP4313Re
doublets_checking(HP4313Re_no_doublets, HP4313Re_with_doublets, "HP4313Re")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_no_doublets.RData")
HP4568Pr_no_doublets <- HP4568Pr
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_with_doublets.RData")
HP4568Pr_with_doublets <- HP4568Pr
doublets_checking(HP4568Pr_no_doublets, HP4568Pr_with_doublets, "HP4568Pr")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_no_doublets.RData")
HP4568Re_no_doublets <- HP4568Re
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_with_doublets.RData")
HP4568Re_with_doublets <- HP4568Re
doublets_checking(HP4568Re_no_doublets, HP4568Re_with_doublets, "HP4568Re")
