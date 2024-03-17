library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(funrar)
library(ggplot2)
library(cowplot)
library(paletteer)
library(rlang)
library(wakefield)
library(vegan)

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_3.R")

### Figure 3 left half ---

## Read in the color legend --
global_genera_cols <- readr::read_rds(here("slccolon/Global_Genera_Cols.RDS"))
names(global_genera_cols)

## Trios --
trios_lc_barplot <- generate_L6_taxa_plots("slccolon/Trios/taxa_barplots/LumCol_level-6_trios.csv","Colon Lumen", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

trios_mc_barplot <- generate_L6_taxa_plots("slccolon/Trios/taxa_barplots/MucCol_level-6.csv","Colon Mucosa", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

## Long Term --
lt_lc_barplot <- generate_L6_taxa_plots("slccolon/Long_Term/taxa_barplots/LumCol_level-6.csv","Colon Lumen", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

lt_mc_barplot <- generate_L6_taxa_plots("slccolon/Long_Term/taxa_barplots/MucCol_level-6.csv","Colon Mucosa", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

## Baseline --
jax_lc_barplot <- generate_L6_taxa_plots("slccolon/Baseline/taxa_barplots/JAX_Baseline_level-6.csv","Fecal Pellet", ".*g__",global_genera_cols) +
  theme(legend.position = "none")

left_half <- plot_grid(trios_lc_barplot, trios_mc_barplot,
                lt_lc_barplot, lt_mc_barplot, 
                jax_lc_barplot, NULL, 
                nrow=3, label_size = 16, labels=c("A","B","C","D","E"))

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

### From archived script Trios_RS_Jensen_Shannon.R
here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_3.R")

metadata <- read.table(here("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv"), header=TRUE)
counts <- read.table(here("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv"), header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon 
trios_lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(trios_lumcol_meta) <- trios_lumcol_meta$SampleID
lumcol <- trios_lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon
trios_muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(trios_muccol_meta) <- trios_muccol_meta$SampleID
muccol <- trios_muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))


## Prevalence filter datasets -- 
# Luminal Colon
trios_lumcol_counts <- prevalence_filter(lumcol_counts,13)

# Mucosal Colon
82*0.15
trios_muccol_counts <- prevalence_filter(muccol_counts,12)


## Calculate RS Jensen Shannon distance matrix -- 
trios_muccol.dist <- calculate_rsjensen(trios_muccol_counts)
trios_lumcol.dist <- calculate_rsjensen(trios_lumcol_counts)

### Baseline ---
metadata <- read.table(here("Baseline/starting_files/Baseline_Metadata.tsv"), header=TRUE)
counts <- read.table(here("Baseline/starting_files/Baseline_ASV_table_Silva_v138_1.tsv"), header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]
JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts))
row.names(JAX_meta) <- JAX_meta$SampleID
f_JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts),Sex=="Female")
m_JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts),Sex=="Male")

JAX_counts <- counts %>% select(all_of(JAX_meta$SampleID))
JAX_counts_prev <- prevalence_filter(JAX_counts,13)
JAX.dist <- calculate_rsjensen(JAX_counts_prev)

### From archived script Long_Term_RS_Jensen_Shannon.R ---
metadata <- read.table("Long_Term/starting_files/SLC_LT_metadata.tsv", header=TRUE)
counts <- read.table("Long_Term/starting_files/SLT_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon 
lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(lumcol_meta) <- lumcol_meta$SampleID
lumcol <- lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon
muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(muccol_meta) <- muccol_meta$SampleID
muccol <- muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))

## Prevalence filter datasets -- 
# Luminal Colon
0.15*90 #13 samples
lumcol_counts <- prevalence_filter(lumcol_counts,13)

# Mucosal Colon 
0.15*89
muccol_counts <- prevalence_filter(muccol_counts,13)

## Calculate RS Jensen Shannon 
muccol.dist <- calculate_rsjensen(muccol_counts)
lumcol.dist <- calculate_rsjensen(lumcol_counts)

### Figure 3 right half ---

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=trios_lumcol.dist,
                                     counts = trios_lumcol_counts,
                                     metadata = trios_lumcol_meta,
                                     title="Colon Lumen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = here("Trios/beta_diversity/LumCol_Top_Taxa_PcoA.csv"))

trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=trios_muccol.dist,
                                     counts = trios_muccol_counts,
                                     metadata = trios_muccol_meta,
                                     title="Colon Mucosa",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = here("Trios/beta_diversity/MucCol_Top_Taxa_PcoA.csv"))

jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=JAX.dist,
                                         counts = JAX_counts_prev,
                                         metadata = JAX_meta,
                                         title="Fecal Pellet",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = here("Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv"))
slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lumcol.dist,
                                   counts = lumcol_counts,
                                   metadata = lumcol_meta,
                                   title="Colon Lumen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = here("Long_Term/LumCol_Top_Taxa_PcoA.csv"))

slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                   counts = muccol_counts,
                                   metadata = muccol_meta,
                                   title="Colon Mucosa",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = here("Long_Term/MucCol_Top_Taxa_PcoA.csv"))
class(slt_mc_pcoa)
class(trios_lc_pcoa)
class(trios_mc_pcoa)
class(slt_lc_pcoa)
class(jax_baseline_pcoa)

#Create right half
right_half <- plot_grid(trios_lc_pcoa, trios_mc_pcoa,
          slt_lc_pcoa, slt_mc_pcoa,
          NULL, jax_baseline_pcoa,nrow = 3,
          labels=c("F","G","H","I","","J"),
          label_size = 16)

dev.new(width=10,height=10)
right_half
#Generate each segment as individual plot
left_half
grid::grid.newpage()
grid::grid.draw(legend)
right_half
