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
                                              "Luminal Colon: MUT vs WT")+
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
                                              "Luminal Colon: MUT vs WT")+
  theme(legend.position = "none")

# Baseline
data<-read.table("Baseline/differential_taxa/JAX_L6_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

jax_mut <- data %>% filter(value=="MUT")


jax_mut_genus <- make_genus_level_taxa_dotplot(jax_mut,
                                               "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                               "Fecal Pellet: MUT vs WT")+
  theme(legend.position = "none")

## Final Figure --
left_half <- plot_grid(tr_genus_mut, lt_genus_mut, jax_mut_genus,
                       nrow=3, 
                       labels=c("A", "B", "C"),
                       label_size=20)

## Mucosal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

trios_l6_het_significant <- data %>% filter(value=="HET")
trios_l6_mut_significant <- data %>% filter(value=="MUT")

tr_genus_mut <- make_genus_level_taxa_dotplot(trios_l6_mut_significant, 
                                              "Trios/differential_taxa/L6_Mucosal_ColonRelative_Abundance-ASV.RDS",
                                              "Trios MucCol MUT: Genus ~ SeqRun + Site + Sex + Genotype + (1|MouseID)")
tr_genus_het <- make_genus_level_taxa_dotplot(trios_l6_het_significant, 
                                              "Trios/differential_taxa/L6_Mucosal_ColonRelative_Abundance-ASV.RDS",
                                              "Trios MucCol HET: Genus ~ SeqRun + Site + Sex + Genotype + (1|MouseID)")

plot_grid(tr_genus_het, tr_genus_mut, labels=c("A", "B"), label_size = 20)

# Long Term 
data<-read.table("Long_Term/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature) 
relA <- readRDS("Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS")
print(row.names(relA))
lt_l6_het_significant <- data %>% filter(value=="HET")
lt_l6_mut_significant <- data %>% filter(value=="MUT")
print(lt_l6_mut_significant$feature)

lt_genus_mut <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                              "Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS",
                                              "LongTerm MucCol MUT: Genus ~ Site + Sex + Genotype + (1|MouseID)")
lt_genus_het <- make_genus_level_taxa_dotplot(lt_l6_het_significant, 
                                              "Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS",
                                              "LongTerm MucCol HET: Genus ~ Site + Sex + Genotype + (1|MouseID)")

plot_grid(lt_genus_het, lt_genus_mut, labels=c("A", "B"), label_size = 20)

