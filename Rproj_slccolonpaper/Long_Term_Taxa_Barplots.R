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

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") 
setwd("/home/julianne/Documents/slccolonpaper/slccolon/")

here::i_am("Rproj_slccolonpaper/Long_Term_Taxa_Barplots.R")

## Replace genera names with legible ones --
wrangle_genera_names("Long_Term/taxa_barplots/LumCol_level-6.csv", "Long_Term/taxa_barplots/","LuminalColon_level-6.RDS")
wrangle_genera_names("Long_Term/taxa_barplots/MucCol_level-6.csv", "Long_Term/taxa_barplots/","MucosalColon_level-6.RDS")

## Plot the barplots --
trios_lc_barplot <- generate_L6_taxa_plots("Trios/LuminalColon_level-6.RDS","Luminal Colon", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

trios_mc_barplot <- generate_L6_taxa_plots("Trios/MucosalColon_level-6.RDS","Mucosal Colon", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

## Draw the legend --
L6_legend <- generate_L6_taxa_plots("Trios/MucosalColon_level-6.RDS","Mucosal Colon", ".*g__",global_genera_cols) +
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
L2_lum<-readRDS("Long_Term/taxa_barplots/LuminalColon_level-6.RDS")
L2_lum<- as.matrix(L2_lum)
L2_lum<-funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
L2_lum <-select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))

L2_lum$Taxa <-taxa
labels_lum <- unique(L2_lum$Taxa)
print(labels_lum)
## Extract taxa >0.001 for use in legend ---
L2_lum<-read.csv("Long_Term/taxa_barplots/MucCol_level-6.csv",header=TRUE,row.names=1)
L2_lum<- as.matrix(L2_lum)
L2_lum<-funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
L2_lum <-select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))
taxa<-gsub(".*g__","",taxa )

L2_lum$Taxa <-taxa
labels_muc <- unique(L2_lum$Taxa)

global <- union(labels_lum,labels_muc)
length(global)
print(global)
## Generate a color key using paletteer colors ---

add_cols2 <- paletteer_d("ggthemes::Classic_20",20)	
add_cols4 <- paletteer_d("ggthemes::calc",12)
add_cols3 <- paletteer_d("ggsci::category20_d3", 20)
add_cols5 <- paletteer_d("basetheme::clean",2)
add_cols6 <- paletteer_d("basetheme::dark",10)
global_genera_cols <- c(add_cols2,add_cols3,add_cols6, add_cols4,add_cols5)
global_genera_cols <- unique(global_genera_cols)
names(global_genera_cols) <- global
saveRDS(global_genera_cols,"Trios/Genera_cols.RDS")




