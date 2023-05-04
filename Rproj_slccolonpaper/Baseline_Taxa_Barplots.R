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

here::i_am("Rproj_slccolonpaper/Baseline_Taxa_Barplots.R")

## Replace genera names with legible ones --
wrangle_genera_names("Baseline/taxa_barplots/JAX_Baseline_level-6.csv", "Baseline/taxa_barplots/","JAX_level-6.RDS")


## Plot the barplots --
jax_lc_barplot <- generate_L6_taxa_plots("Baseline/taxa_barplots/JAX_Baseline_level-6.csv","Fecal Pellet", ".*g__",global_genera_cols) +
  theme(legend.position = "none")



## Draw the legend --
L6_legend <- generate_L6_taxa_plots("Baseline/taxa_barplots/JAX_Baseline_level-6.csv","Fecal Pellets", ".*g__",global_genera_cols) +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_legend)
grid::grid.newpage()
grid::grid.draw(legend)

## Final Figures --
plot_grid(trios_lc_barplot,trios_mc_barplot)


## Extract taxa >0.01 for use in legend ---
L2_lum<-read.csv("Baseline/taxa_barplots/JAX_Baseline_level-6.csv",row.names=1, header=TRUE)
L2_lum<- as.matrix(L2_lum)
L2_lum<-funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.01, row.names(L2_lum), "Other"))
L2_lum <-select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
print(taxa)
L2_lum <- select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))
family<-gsub(".*f__","",taxa )
genus <- gsub(".*g__","",taxa)
order <- gsub(".*o__","",taxa)

L2_lum$Family<-gsub(".g__.*","",family)
L2_lum$Genus <- genus
L2_lum$Order<-gsub(".f__.*","",order)

L2_lum <- L2_lum %>% mutate(annotation = ifelse(L2_lum$Genus=="", paste0(L2_lum$Family,"..f."), L2_lum$Genus))
L2_lum <- L2_lum %>% mutate(annotation = ifelse(L2_lum$Family=="", paste0(L2_lum$Order,"..o."), L2_lum$annotation))

labels_muc <- unique(L2_lum$annotation)
print(taxa)
jax_global <- union(labels_lum,labels_muc)
length(jax_global)

## Generate a color key using paletteer colors ---

lt_trios_jax_global <- c(lt_global, trios_global,jax_global)
lt_trios_jax_global <- unique(lt_trios_jax_global)
add_cols2 <- paletteer_d("ggthemes::Classic_20",20)	
add_cols4 <- paletteer_d("ggthemes::calc",8)
global_genera_cols <- c(add_cols2, add_cols4) 
global_genera_cols <- unique(global_genera_cols)

saveRDS(global_genera_cols,"Global_Genera_Cols.RDS")




