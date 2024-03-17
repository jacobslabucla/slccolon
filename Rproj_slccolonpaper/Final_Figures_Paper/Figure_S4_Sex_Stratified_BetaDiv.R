library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(funrar)
library(ggplot2)
library(cowplot)
library(paletteer)
library(rlang)
library(wakefield)
library(vegan)

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_S4_Sex_Stratified_BetaDiv.R")

### From archived script Trios_RS_Jensen_Shannon.R
here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_3.R")

metadata <- read.table(here("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv"), header=TRUE)
counts <- read.table(here("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv"), header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

# Luminal Colon F
f_trios_lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", Sex=="Female",SampleID %in% names(counts))
row.names(f_trios_lumcol_meta) <- f_trios_lumcol_meta$SampleID
f_lumcol <- f_trios_lumcol_meta$SampleID
f_lumcol_counts <- counts %>% select(all_of(f_lumcol))

# Luminal Colon M
m_trios_lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", Sex=="Male",SampleID %in% names(counts))
row.names(m_trios_lumcol_meta) <- m_trios_lumcol_meta$SampleID
m_lumcol <- m_trios_lumcol_meta$SampleID
m_lumcol_counts <- counts %>% select(all_of(m_lumcol))

# Mucosal Colon F
f_trios_muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", Sex=="Female",SampleID %in% names(counts))
row.names(f_trios_muccol_meta) <- f_trios_muccol_meta$SampleID
f_muccol <- f_trios_muccol_meta$SampleID
f_muccol_counts <- counts %>% select(all_of(f_muccol))

# Mucosal Colon M
m_trios_muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", Sex=="Male",SampleID %in% names(counts))
row.names(m_trios_muccol_meta) <- m_trios_muccol_meta$SampleID
m_muccol <- m_trios_muccol_meta$SampleID
m_muccol_counts <- counts %>% select(all_of(m_muccol))

# Luminal Colon F
36*0.15
f_trios_lumcol_counts <- prevalence_filter(f_lumcol_counts,5)

# Luminal Colon M
53*0.15
m_trios_lumcol_counts <- prevalence_filter(m_lumcol_counts,8)

# Mucosal Colon M
51*0.15
m_trios_muccol_counts <- prevalence_filter(m_muccol_counts,8)

# Mucosal Colon F
53*0.15
f_trios_muccol_counts <- prevalence_filter(f_muccol_counts,8)

## Calculate RS Jensen Shannon distance matrix -- 
f_trios_muccol.dist <- calculate_rsjensen(f_trios_muccol_counts)
m_trios_muccol.dist <- calculate_rsjensen(m_trios_muccol_counts)
f_trios_lumcol.dist <- calculate_rsjensen(f_trios_lumcol_counts)
m_trios_lumcol.dist <- calculate_rsjensen(m_trios_lumcol_counts)

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

f_JAX_counts <- counts %>% select(all_of(f_JAX_meta$SampleID))
f_JAX_counts_prev <- prevalence_filter(f_JAX_counts,6)
f_JAX.dist <- calculate_rsjensen(f_JAX_counts_prev)

m_JAX_counts <- counts %>% select(all_of(m_JAX_meta$SampleID))
m_JAX_counts_prev <- prevalence_filter(m_JAX_counts,7)
m_JAX.dist <- calculate_rsjensen(m_JAX_counts_prev)

### Figure S4 ---
f_trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=f_trios_lumcol.dist,
                                     counts = f_trios_lumcol_counts,
                                     metadata = f_trios_lumcol_meta,
                                     title="2- Month Females: Colon Lumen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = here("Trios/beta_diversity/F_LumCol_Top_Taxa_PcoA.csv"))

m_trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=m_trios_lumcol.dist,
                                       counts = m_trios_lumcol_counts,
                                       metadata = m_trios_lumcol_meta,
                                       title="2- Month Males: Colon Lumen",
                                       colorvariable = Genotype,
                                       colorvector = cols,
                                       wa_scores_filepath = here("Trios/beta_diversity/M_LumCol_Top_Taxa_PcoA.csv"))

f_trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=f_trios_muccol.dist,
                                     counts = f_trios_muccol_counts,
                                     metadata = f_trios_muccol_meta,
                                     title="2- Month Females: Colon Mucosa",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = here("Trios/beta_diversity/f_MucCol_Top_Taxa_PcoA.csv"))

m_trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=m_trios_muccol.dist,
                                       counts = m_trios_muccol_counts,
                                       metadata = m_trios_muccol_meta,
                                       title="2- Month Males: Colon Mucosa",
                                       colorvariable = Genotype,
                                       colorvector = cols,
                                       wa_scores_filepath = here("Trios/beta_diversity/m_MucCol_Top_Taxa_PcoA.csv"))

f_jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=f_JAX.dist,
                                         counts = f_JAX_counts_prev,
                                         metadata = f_JAX_meta,
                                         title="3-4 month Females: Fecal Pellet",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = here("Baseline/beta_diversity/f_JAX_Top_Taxa_PcoA.csv"))

m_jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=m_JAX.dist,
                                           counts = m_JAX_counts_prev,
                                           metadata = m_JAX_meta,
                                           title="3-4 month Males: Fecal Pellet",
                                           colorvariable = Genotype,
                                           colorvector = cols,
                                           wa_scores_filepath = here("Baseline/beta_diversity/m_JAX_Top_Taxa_PcoA.csv"))

dev.new(width=10,height=10)
plot_grid(f_trios_mc_pcoa, m_trios_mc_pcoa,
          f_trios_lc_pcoa, m_trios_lc_pcoa,
          f_jax_baseline_pcoa,m_jax_baseline_pcoa,
          labels=c("A","B","C","D","E","F"),nrow=3)

### Statistics ---

## Trios --
# Luminal Colon Females
data.dist<-f_trios_lumcol.dist
metadata <- f_trios_lumcol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Litter +  Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Luminal Colon Males
data.dist<-m_trios_lumcol.dist
metadata <- m_trios_lumcol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Litter +  Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon Females
data.dist<-f_trios_muccol.dist
metadata <- f_trios_muccol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Litter + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon Males
data.dist<-m_trios_muccol.dist
metadata <- m_trios_muccol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Litter + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

## Baseline --
# Females
data.dist<-f_JAX.dist
metadata <- f_JAX_meta
row.names(metadata) <- f_JAX_meta$SampleID

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Males
data.dist<-m_JAX.dist
metadata <- m_JAX_meta
row.names(metadata) <- m_JAX_meta$SampleID

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab