library(ggplot2)
library(cowplot)
library(dplyr)

here::i_am("Rproj_slccolonpaper/Overlapping_Genera.R")

## Luminal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Luminal_SI_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/all_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

trios_l6_het_significant <- data %>% filter(value=="HET")
trios_l6_mut_significant <- data %>% filter(value=="MUT")

tr_genus_mut <- make_genus_level_taxa_dotplot(trios_l6_mut_significant, 
                                              "Trios/differential_taxa/L6_Luminal_ColonRelative_Abundance-ASV.RDS",
                                              "Trios LumCol MUT: Genus ~ SeqRun + Site + Sex + Genotype + (1|MouseID)")
tr_genus_het <- make_genus_level_taxa_dotplot(trios_l6_het_significant, 
                                              "Trios/differential_taxa/L6_Luminal_ColonRelative_Abundance-ASV.RDS",
                                              "Trios LumCol HET: Genus ~ SeqRun + Site + Sex + Genotype + (1|MouseID)")

plot_grid(tr_genus_het, tr_genus_mut, labels=c("A", "B"), label_size = 20)

# Long Term 
data<-read.table("Long_Term/differential_taxa/L6_Luminal_SI_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/all_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)
relA <- readRDS("Long_Term/differential_taxa/L6_Luminal_SI-ASV.RDS")
print(row.names(relA))
lt_l6_het_significant <- data %>% filter(value=="HET")
lt_l6_mut_significant <- data %>% filter(value=="MUT")
print(lt_l6_mut_significant$feature)

lt_genus_mut <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                              "Long_Term/differential_taxa/L6_Luminal_SI-ASV.RDS",
                                              "Luminal SI")
lt_genus_het <- make_genus_level_taxa_dotplot(lt_l6_het_significant, 
                                              "Long_Term/differential_taxa/L6_Luminal_SI-ASV.RDS",
                                              "LongTerm LumSI HET: Genus ~ Site + Sex + Genotype + (1|MouseID)")

plot_grid(lt_genus_het, lt_genus_mut, labels=c("A", "B"), label_size = 20)

# Baseline
data<-read.table("Baseline/differential_taxa/JAX_L6_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

jax_mut <- data %>% filter(value=="MUT")
jax_het <- data %>% filter(value=="HET")


jax_mut_genus <- make_genus_level_taxa_dotplot(jax_mut,
                                               "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                               "Baseline MUT: Genus~ Sex + Genotype")


jax_het_genus <- make_genus_level_taxa_dotplot(jax_het,
                                               "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                               "Baseline HET: Genus~ Sex + Genotype")

plot_grid(jax_het_genus, jax_mut_genus, labels= c("A","B"), label_size = 20)

## Mucosal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Mucosal_SI_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/all_results.tsv", header=TRUE)
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
data<-read.table("Long_Term/differential_taxa/L6_Mucosal_SI_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature) 
relA <- readRDS("Long_Term/differential_taxa/L6_Mucosal_SI-ASV.RDS")
print(row.names(relA))
lt_l6_het_significant <- data %>% filter(value=="HET")
lt_l6_mut_significant <- data %>% filter(value=="MUT")
print(lt_l6_mut_significant$feature)

lt_genus_mut <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                              "Long_Term/differential_taxa/L6_Mucosal_SI-ASV.RDS",
                                              "Mucosal SI")
lt_genus_het <- make_genus_level_taxa_dotplot(lt_l6_het_significant, 
                                              "Long_Term/differential_taxa/L6_Mucosal_SI-ASV.RDS",
                                              "LongTerm MucSI HET: Genus ~ Site + Sex + Genotype + (1|MouseID)")

plot_grid(lt_genus_het, lt_genus_mut, labels=c("A", "B"), label_size = 20)

