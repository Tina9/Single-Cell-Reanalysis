library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(ggplot2)
library(SingleR)
options(scipen = 20)

ClusterFinding <- function(pbmc, sample_name){
  mouse_refer <- ImmGenData()
  express_data <- GetAssayData(pbmc, slot = "data")
  metadata <- pbmc@meta.data
  annote_res <- SingleR(test = express_data,
                        ref = mouse_refer,
                        labels = mouse_refer$label.main)
  metadata$Immlabels <- as.data.frame(annote_res)$labels
  cluster_identity_file <- paste(sample_name, "_markers_no_doublets_ref.txt", sep = "")
  write.table(metadata, file = cluster_identity_file, sep = "\t", row.names = F)
  
  pbmc@meta.data <- metadata
  return(pbmc)
}

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
HP4540Pr <- ClusterFinding(HP4540Pr, "HP4540Pr")
save(HP4540Pr, file = "HP4540Pr_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_no_doublets.RData")
HP4313Pr <- ClusterFinding(HP4313Pr, "HP4313Pr")
save(HP4313Pr, file = "HP4313Pr_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_no_doublets.RData")
HP4313Re <- ClusterFinding(HP4313Re, "HP4313Re")
save(HP4313Re, file = "HP4313Re_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_no_doublets.RData")
HP4568Pr <- ClusterFinding(HP4568Pr, "HP4568Pr")
save(HP4568Pr, file = "HP4568Pr_no_doublets.RData")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_no_doublets.RData")
HP4568Re <- ClusterFinding(HP4568Re, "HP4568Re")
save(HP4568Re, file = "HP4568Re_no_doublets.RData")