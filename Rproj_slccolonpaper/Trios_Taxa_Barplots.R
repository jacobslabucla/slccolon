### Trios - Taxa summary plots 
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(funrar)
library(ggplot2)
library(cowplot)
library(paletteer)
library(rlang)
library(tidyverse)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
setwd("../") # CHANGE to the directory containing the fastq files

here::i_am("Rproj_slccolonpaper/Trios_Taxa_Barplots.R")

metadata <- read.table("Trios/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
lumcol_counts <- read.csv("Trios/LumCol_level-6_trios.csv", header = TRUE, row.names=1)

muccol_counts <- read.csv("Trios/MucCol_level-6.csv", header = TRUE, row.names=1)

## Replace genera names with legible ones --
wrangle_genera_names("Trios/LumCol_level-6_trios.csv", "Trios/","LuminalColon_level-6.RDS")
wrangle_genera_names("Trios/MucCol_level-6.csv", "Trios/","MucosalColon_level-6.RDS")

## Plot the barplots --
trios_lc_barplot <- generate_L6_taxa_plots("Trios/taxa_barplots/LuminalColon_level-6.RDS","Luminal Colon", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

trios_mc_barplot <- generate_L6_taxa_plots("Trios/taxa_barplots/MucosalColon_level-6.RDS","Mucosal Colon", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

plot_grid(trios_lc_barplot, trios_mc_barplot)

## Draw the legend --
L6_legend <- generate_L6_taxa_plots("Trios/taxa_barplots/MucosalColon_level-6.RDS","Mucosal Colon", ".*g__",global_genera_cols) +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_legend)
grid::grid.newpage()
grid::grid.draw(legend)

## Final Figures --
plot_grid(trios_lc_barplot,trios_mc_barplot)
## Extract taxa >0.001 for use in legend ---
L2_lum<-read_rds("Trios/taxa_barplots/LuminalColon_level-6.RDS")
L2_lum<- as.matrix(L2_lum)
L2_lum<-funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.01, row.names(L2_lum), "Other"))
L2_lum <-select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))


L2_lum$Taxa <-taxa
labels_lum <- unique(L2_lum$Taxa)

## Extract taxa >0.001 for use in legend ---
L2_lum<-read.csv("Trios/taxa_barplots/MucCol_level-6.csv",header=TRUE,row.names=1)
L2_lum<- as.matrix(L2_lum)
L2_lum<-funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.01, row.names(L2_lum), "Other"))
L2_lum <-select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))
taxa<-gsub(".*g__","",taxa )

L2_lum$Taxa <-taxa
labels_muc <- unique(L2_lum$Taxa)

trios_global <- union(labels_lum,labels_muc)
length(global)
## Generate a color key using paletteer colors ---

trios_cols <- paletteer_d("ggthemes::Classic_20",20)	

global_genera_cols <- unique(trios_cols)
names(global_genera_cols) <- global
saveRDS(global_genera_cols,"Trios/Genera_cols.RDS")




