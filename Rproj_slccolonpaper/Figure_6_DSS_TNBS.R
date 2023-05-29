library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)

setwd("/home/julianne/Documents/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Figure_6_DSS_TNBS.R")

plot_grid(percent_weight, NULL, dss_histo_plot,
          percent_weight_tnbs, NULL, colon_plot,
          nrow=2, label_size = 20,
          labels=c("A", "B","C","D","E","F"))
