library(ggplot2)
library(cowplot)
library(dplyr)
library(paletteer)

# rm(list = ls())
here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_4.R")

## generate phylum color vector 
phyla_names <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Deinococcota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", 
                 "Deferribacterota", "Cyanobacteria")
phyla_colors <- paletteer_d("ggthemes::calc",9)
names(phyla_colors) <- phyla_names

## Luminal --
# Trios - Female Lum Col
data<-readr::read_delim(here("Trios/differential_taxa/F_L6_Luminal_Colon_L6_Maaslin2_Site_Genotype_1-MouseID/significant_results.tsv"))
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

f_trios_l6_mut_significant <- data %>% filter(value=="MUT")

f_tr_genus_mut <- make_genus_level_taxa_dotplot(f_trios_l6_mut_significant, 
                                              "Trios/differential_taxa/F_L6_Luminal_ColonRelative_Abundance-ASV.RDS",
                                              "Female Colon Lumen",
                                              phyla_colors)+
  theme(legend.position = "none")

# Trios - Male Lum Col
data<-readr::read_delim(here("Trios/differential_taxa/M_L6_Luminal_Colon_L6_Maaslin2_Sequencing_Run_Site_Genotype_1-MouseID/significant_results.tsv"))
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

m_trios_l6_mut_significant <- data %>% filter(value=="MUT")

m_tr_genus_mut <- make_genus_level_taxa_dotplot(m_trios_l6_mut_significant, 
                                                "Trios/differential_taxa/M_L6_Luminal_ColonRelative_Abundance-ASV.RDS",
                                                "Male Colon Lumen",
                                                phyla_colors)+
  theme(legend.position = "none")

# Trios- Mucosal Colon F
data<-readr::read_delim(here("Trios/differential_taxa/F_L6_Mucosal_Colon_L6_Maaslin2_Site_Genotype_1-MouseID/significant_results.tsv"))
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

f_tr_mc_l6_mut_significant <- data %>% filter(value=="MUT")
print(f_tr_mc_l6_mut_significant$feature)

f_tr_mc_genus_mut <- make_genus_level_taxa_dotplot(f_tr_mc_l6_mut_significant, 
                                              "Trios/differential_taxa/F_L6_Mucosal_ColonRelative_Abundance-ASV.RDS",
                                              "Female Colon Mucosa",phyla_colors)+
  theme(legend.position = "none")

# Trios- Mucosal Colon M
data<-readr::read_delim(here("Trios/differential_taxa/M_L6_Mucosal_Colon_L6_Maaslin2_Sequencing_Run_Site_Genotype_1-MouseID/significant_results.tsv"))
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

m_tr_mc_l6_mut_significant <- data %>% filter(value=="MUT")
print(m_tr_mc_l6_mut_significant$feature)

m_tr_mc_genus_mut <- make_genus_level_taxa_dotplot(m_tr_mc_l6_mut_significant, 
                                                   "Trios/differential_taxa/M_L6_Mucosal_ColonRelative_Abundance-ASV.RDS",
                                                   "Male Colon Mucosa",phyla_colors)+
  theme(legend.position = "none")

# Baseline F
data<-readr::read_delim(here("Baseline/differential_taxa/f_JAX_L6_Maaslin2_Genotype/significant_results.tsv"))
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

f_jax_mut <- data %>% filter(value=="MUT")
print(f_jax_mut)

f_jax_mut_genus <- make_genus_level_taxa_dotplot(f_jax_mut,
                                               "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                               "Female Fecal Pellet",phyla_colors)+
  theme(legend.position = "none")

# Baseline M
data<-readr::read_delim(here("Baseline/differential_taxa/m_JAX_L6_Maaslin2_Genotype/significant_results.tsv"))
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

m_jax_mut <- data %>% filter(value=="MUT")
print(m_jax_mut)

m_jax_mut_genus <- make_genus_level_taxa_dotplot(m_jax_mut,
                                                 "Baseline/differential_taxa/JAX_L6Relative_Abundance_ASV.RDS",
                                                 "Male Fecal Pellet",phyla_colors)+
  theme(legend.position = "none")

## Get legend 
L2_legend <- make_genus_level_taxa_dotplot(lt_l6_mut_significant, 
                                           "Long_Term/differential_taxa/L6_Mucosal_Colon-ASV.RDS",
                                           "Colon Mucosa", phyla_colors)+
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_legend)

## Final Plots --

left_half <- plot_grid(f_tr_genus_mut, m_tr_genus_mut, f_jax_mut_genus,
                       nrow=3, 
                       labels=c("A", "B", "C"),
                       label_size=20,
                       rel_heights = c(0.5,1, 0.5))
dev.new(width=10,height=10)
left_half

right_half <- plot_grid(f_tr_mc_genus_mut, m_tr_mc_genus_mut,
                        m_jax_mut_genus,
                        nrow=3, 
                        labels=c("D", "E","F"),
                        label_size=20,
                        rel_heights = c(1,1,1))
dev.new(width=10,height=10)
right_half 

# Arranged plots
plot_grid(plotlist = list(left_half, right_half))

# Legend
grid::grid.newpage()
grid::grid.draw(legend)
