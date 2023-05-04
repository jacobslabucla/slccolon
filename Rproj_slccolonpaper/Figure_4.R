library(ggplot2)
library(cowplot)
library(dplyr)

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
                                              "Luminal Colon")+
  theme(legend.position = "none")




# Long Term 
data<-read.table("Long_Term/differential_taxa/L6_Luminal_Colon_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

lt_l6_mut_significant <- data %>% filter(value=="MUT")
print(lt_l6_mut_significant$feature)

lt_genus_mut <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                              "Long_Term/differential_taxa/L6_Luminal_Colon-ASV.RDS",
                                              "Luminal Colon")+
  theme(legend.position = "none")

# Baseline
data<-read.table("Baseline/differential_taxa/JAX_L6_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

jax_mut <- data %>% filter(value=="MUT")


jax_mut_genus <- make_genus_level_taxa_dotplot(jax_mut,
                                               "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                               "Fecal Pellet")+
  theme(legend.position = "none")



## Mucosal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)


trios_l6_mut_significant <- data %>% filter(value=="MUT")

muc_tr_genus_mut <- make_genus_level_taxa_dotplot(trios_l6_mut_significant, 
                                              "Trios/differential_taxa/L6_Mucosal_ColonRelative_Abundance-ASV.RDS",
                                              "Colon Mucosa") + 
                    theme(legend.position = "none")

# Long Term 
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
                                              "Colon Mucosa")+
                    theme(legend.position="none")

## Final Figure --
left_half <- plot_grid(tr_genus_mut, lt_genus_mut, jax_mut_genus,
                       nrow=3, 
                       labels=c("A", "B", "C"),
                       label_size=20)

right_half <- plot_grid(muc_tr_genus_mut, muc_lt_genus_mut,
                        nrow=2, 
                        labels=c("C", "D"),
                        label_size=20,
                        rel_heights = c(0.5,1))
