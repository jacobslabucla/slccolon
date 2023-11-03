library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(funrar)
library(ggplot2)
library(cowplot)
library(paletteer)
library(rlang)
library(wakefield)
library(vegan)

here::i_am("Rproj_slccolonpaper/Figure_3.R")

### Figure 3 left half ---

## Read in the color legend --
global_genera_cols <- readr::read_rds(here("Global_Genera_Cols.RDS"))
names(global_genera_cols)

## Trios --
trios_lc_barplot <- generate_L6_taxa_plots("Trios/taxa_barplots/LumCol_level-6_trios.csv","Colon Lumen", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

trios_mc_barplot <- generate_L6_taxa_plots("Trios/taxa_barplots/MucCol_level-6.csv","Colon Mucosa", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

plot_grid(trios_lc_barplot, trios_mc_barplot)


## Long Term --
lt_lc_barplot <- generate_L6_taxa_plots("Long_Term/taxa_barplots/LumCol_level-6.csv","Colon Lumen", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

lt_mc_barplot <- generate_L6_taxa_plots("Long_Term/taxa_barplots/MucCol_level-6.csv","Colon Mucosa", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

plot_grid(lt_lc_barplot, lt_mc_barplot)

## Baseline --
jax_lc_barplot <- generate_L6_taxa_plots("Baseline/taxa_barplots/JAX_Baseline_level-6.csv","Fecal Pellet", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

left_half <- plot_grid(trios_lc_barplot, trios_mc_barplot,
                lt_lc_barplot, lt_mc_barplot, 
                jax_lc_barplot, NULL, 
                nrow=3, label_size = 16, labels=c("A","B","C","D","E"))

left_half

## Legend --
## Draw the legend --
fake_df <- data.frame(taxa = names(global_genera_cols),
                      var1 = seq(1,28,by=1),
                      var2 = upper(n=28, k = 5, x = LETTERS, prob = NULL, name = "Upper"))
L6_legend <- ggplot(fake_df, aes(x=var2, y=var1, color=taxa,fill=taxa)) +
  geom_point(pch=15,size=4,position=position_jitter(width=0.25),alpha=1, aes(color=taxa, fill=taxa))+
  scale_color_manual(values=global_genera_cols) +
  scale_fill_manual(values=global_genera_cols) +
  theme_cowplot(16)+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=28, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_legend)
grid::grid.newpage()
grid::grid.draw(legend)
                                                                                                                                              

### Figure 3 right half ---

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=trios_lumcol.dist,
                                     counts = trios_lumcol_counts,
                                     metadata = trios_lumcol_meta,
                                     title="Colon Lumen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/LumCol_Top_Taxa_PcoA.csv")

trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=trios_muccol.dist,
                                     counts = trios_muccol_counts,
                                     metadata = trios_muccol_meta,
                                     title="Colon Mucosa",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/MucCol_Top_Taxa_PcoA.csv")

jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=JAX.dist,
                                         counts = JAX_counts_prev,
                                         metadata = JAX_meta,
                                         title="Fecal Pellet",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = "Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv")
slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lumcol.dist,
                                   counts = lumcol_counts,
                                   metadata = lumcol_meta,
                                   title="Colon Lumen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/LumCol_Top_Taxa_PcoA.csv")

slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                   counts = muccol_counts,
                                   metadata = muccol_meta,
                                   title="Colon Mucosa",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/MucCol_Top_Taxa_PcoA.csv")
class(slt_mc_pcoa)
class(trios_lc_pcoa)
class(trios_mc_pcoa)
class(slt_lc_pcoa)
class(jax_baseline_pcoa)

right_half <- plot_grid(trios_lc_pcoa, trios_mc_pcoa,
          slt_lc_pcoa, slt_mc_pcoa,
          NULL, jax_baseline_pcoa,nrow = 3,
          labels=c("F","G","H","I","","J"),
          label_size=16)
