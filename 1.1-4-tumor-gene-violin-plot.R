library(reshape2)
library(ggplot2)

chosen_gene_substract <- function(pbmc, chosen_genes, sample_name){
  rownames(pbmc@meta.data) <- colnames(GetAssayData(pbmc))
  pbmc@meta.data$identified_clusters[pbmc@meta.data$identified_clusters == "Tumor-proliferating"] <- "Tumor"
  Tumor_data <- subset(pbmc, subset = nFeature_RNA > 500 &
                         identified_clusters == "Tumor")
  scaled_data <- Tumor_data[["RNA"]]@scale.data
  exist_genes <- rownames(scaled_data)[rownames(scaled_data) %in% chosen_genes]
  Tumor_chosen_data <- t(scaled_data[exist_genes,])
  plot_data <- melt(Tumor_chosen_data)
  colnames(plot_data) <- c("barcode", "Gene", "ExpValue")
  plot_data$sample <- sample_name
  return(plot_data)
}

for(i in chosen_genes){
  chosen_gene <- Plot_data[Plot_data$Gene == i,]
  p1 <- ggplot(chosen_gene, aes(x = sample, y = ExpValue, fill = sample, colour = sample)) +
    geom_violin() +
    #geom_jitter(shape = 16, position = position_jitter(0.2)) +
    theme_classic()
  pdf_name <- paste(i, "_exp_vln.pdf", sep = "")
  pdf(pdf_name, height = 9, width = 9)
  plot(p1)
  dev.off()
}

chosen_genes <- c("Mki67", "Birc5", "Cdk1", "Ccnd1", "Pcna", 
                  "Krt14", "Krt17", "Krt8", "Krt18", "Cdh1", 
                  "Epcam", "Cd24a","Vim", "Fn1", "Cdh2")

#setwd("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr")
load("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
HP4540Pr_Tumor_Plot_data <- chosen_gene_substract(HP4540Pr, chosen_genes, "HP4540Pr")

load("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Re/HP4540Re2_no_doublets.RData")
HP4540Re_Tumor_Plot_data <- chosen_gene_substract(HP4540Re2, chosen_genes, "HP4540Re")

load("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_no_doublets.RData")
HP4313Pr_Tumor_Plot_data <- chosen_gene_substract(HP4313Pr, chosen_genes, "HP4313Pr")

load("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_no_doublets.RData")
HP4313Re_Tumor_Plot_data <- chosen_gene_substract(HP4313Re, chosen_genes, "HP4313Re")

load("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_no_doublets.RData")
HP4568Pr_Tumor_Plot_data <- chosen_gene_substract(HP4568Pr, chosen_genes, "HP4568Pr")

load("/Users/xuzhang/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_no_doublets.RData")
HP4568Re_Tumor_Plot_data <- chosen_gene_substract(HP4568Re, chosen_genes, "HP4568Re")

Plot_data <- rbind(HP4540Pr_Tumor_Plot_data,
                   HP4540Re_Tumor_Plot_data,
                   HP4313Pr_Tumor_Plot_data,
                   HP4313Re_Tumor_Plot_data,
                   HP4568Pr_Tumor_Plot_data,
                   HP4568Re_Tumor_Plot_data)

