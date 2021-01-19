library(dplyr)
library(RColorBrewer)

component_stat <- function(pbmc, sample_name){
  #pbmc <- HP4568Pr
  metadata <- pbmc@meta.data
  plot_stat <- metadata %>% 
    dplyr::count(identified_clusters) %>%
    as.data.frame()
  plot_percent_stat <- plot_stat %>%
    mutate(per = (n/sum(n) * 100))
  plot_percent_stat$orig.ident <- sample_name
  return(plot_percent_stat)
}

percent_plot <- function(pbmc_plot, sample_name){
  
  nb.cols <- length(unique(pbmc_plot$identified_clusters))
  mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

  p1 <- ggplot(data = pbmc_plot, aes(x = orig.ident, y = per, fill = identified_clusters)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = mycolors) +
    #facet_grid(~sample) +
    theme_classic() +
    xlab("Samples") +
    ylab("Cluster Percent") +
    labs(fill = "Group", size = 12) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(colour = "black")
    )
  
  barplot_file <- paste(sample_name, "_percent_barplot.pdf", sep = "")
  pdf(barplot_file, height = 9, width = 9)
  plot(p1)
  dev.off()
}

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Pr/HP4568Pr_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4568Re/HP4568Re_no_doublets.RData")
HP4568Pr_percent_plot <- component_stat(HP4568Pr, "HP4568Pr")
HP4568Re_percent_plot <- component_stat(HP4568Re, "HP4568Re")
HP4568_percent_plot <- rbind(HP4568Pr_percent_plot, HP4568Re_percent_plot)
percent_plot(HP4568_percent_plot, "HP4568")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Pr/HP4313Pr_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4313Re/HP4313Re_no_doublets.RData")
HP4313Pr_percent_plot <- component_stat(HP4313Pr, "HP4313Pr")
HP4313Re_percent_plot <- component_stat(HP4313Re, "HP4313Re")
HP4313_percent_plot <- rbind(HP4313Pr_percent_plot, HP4313Re_percent_plot)
percent_plot(HP4313_percent_plot, "HP4313")

setwd("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Pr/HP4540Pr_no_doublets.RData")
load("~/Dropbox/BRCA1-PARPi-10X/Latest-Version/HP4540Re/HP4540Re2_no_doublets.RData")
HP4540Pr_percent_plot <- component_stat(HP4540Pr, "HP4540Pr")
HP4540Re_percent_plot <- component_stat(HP4540Re2, "HP4540Re")
HP4540_percent_plot <- rbind(HP4540Pr_percent_plot, HP4540Re_percent_plot)
percent_plot(HP4540_percent_plot, "HP4540")
