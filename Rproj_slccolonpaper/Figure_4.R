library(ggplot2)
library(cowplot)
library(dplyr)
library(paletteer)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolon/")
here::i_am("Rproj_slccolonpaper/Figure_4.R")

## Luminal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Luminal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

trios_l6_mut_significant <- data %>% filter(value=="MUT")

tr_genus_mut <- make_genus_level_taxa_dotplot(trios_l6_mut_significant, 
                                              "Trios/differential_taxa/L6_Luminal_ColonRelative_Abundance-ASV.RDS",
                                              "Luminal Colon",
                                              phyla_colors)+
  theme(legend.position = "none")

tr_genus_mut + theme(legend.position="right")


# Long Term 
data<-read.table("Long_Term/differential_taxa/L6_Luminal_Colon_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

lt_l6_mut_significant <- data %>% filter(value=="MUT")
print(lt_l6_mut_significant$feature)

lt_genus_mut <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                              "Long_Term/differential_taxa/L6_Luminal_Colon-ASV.RDS",
                                              "Luminal Colon",phyla_colors)+
  theme(legend.position = "none")
lt_genus_mut + theme(legend.position = "right")

# Baseline
data<-read.table("Baseline/differential_taxa/JAX_L6_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

jax_mut <- data %>% filter(value=="MUT")


jax_mut_genus <- make_genus_level_taxa_dotplot(jax_mut,
                                               "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                               "Fecal Pellet",phyla_colors)+
  theme(legend.position = "none")

jax_mut_genus + theme(legend.position="right")

## Mucosal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)


trios_l6_mut_significant <- data %>% filter(value=="MUT")

muc_tr_genus_mut <- make_genus_level_taxa_dotplot(trios_l6_mut_significant, 
                                              "Trios/differential_taxa/L6_Mucosal_ColonRelative_Abundance-ASV.RDS",
                                              "Colon Mucosa", phyla_colors) + 
                    theme(legend.position = "none") 
  
muc_tr_genus_mut + theme(legend.position = "right")
# Long Term - Fig 4D
data<-read.table("Long_Term/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature) 
relA <- readRDS("Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS")
print(row.names(relA))

lt_l6_mut_significant <- data %>% filter(value=="MUT")
print(lt_l6_mut_significant$feature)

muc_lt_genus_mut <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                              "Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS",
                                              "Colon Mucosa", phyla_colors)+
                    theme(legend.position="none")
muc_lt_genus_mut + theme(legend.position="bottom") guides(fill=guide_legend(nrow=22, byrow=TRUE))+

## Get legend 
L2_legend <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                           "Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS",
                                           "Colon Mucosa", phyla_colors)+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid::grid.newpage()
grid::grid.draw(legend)

## Final Figure --
left_half <- plot_grid(tr_genus_mut, lt_genus_mut, jax_mut_genus,
                       nrow=3, 
                       labels=c("A", "B", "C"),
                       label_size=20,
                       rel_heights = c(0.5,1, 0.5))

right_half <- plot_grid(muc_tr_genus_mut, muc_lt_genus_mut,
                        nrow=2, 
                        labels=c("C", "D"),
                        label_size=20,
                        rel_heights = c(0.5,1))

## generate phylum color vector 
phyla_names <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Deinococcota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", 
                 "Deferribacterota", "Cyanobacteria")
phyla_colors <- paletteer_d("ggthemes::calc",9)
names(phyla_colors) <- phyla_names
