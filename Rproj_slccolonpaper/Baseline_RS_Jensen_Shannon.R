library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Baseline_RS_Jensen_Shannon.R")

metadata <- read.table("Baseline/starting_files/Baseline_Metadata.tsv", header=TRUE)
counts <- read.table("Baseline/starting_files/Baseline_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# JAX mice 
JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts))
row.names(JAX_meta) <- JAX_meta$SampleID
JAX <- JAX_meta$SampleID
JAX_counts <- counts %>% select(all_of(JAX))

nohet_JAX_meta <- metadata %>%
  filter(Genotype!="HET") %>%
  filter(Background=="JAX", SampleID %in% names(counts))
row.names(nohet_JAX_meta) <- nohet_JAX_meta$SampleID
JAX <- nohet_JAX_meta$SampleID
nohet_JAX_counts <- counts %>% select(all_of(JAX))

# full dataset
baseline_meta <- metadata
row.names(baseline_meta)=baseline_meta$SampleID
baseline_counts <- counts 

# NIH mice 
NIH_meta <- metadata %>% filter(Background=="NIH", SampleID %in% names(counts))
row.names(JAX_meta) <- JAX_meta$SampleID
NIH <- NIH_meta$SampleID
NIH_counts <- counts %>% select(all_of(NIH))

## Prevalence filter datasets -- 
# JAX
0.15*87
JAX_counts_prev <- prevalence_filter(JAX_counts,13)

0.15*57
nohet_JAX_counts_prev <- prevalence_filter(nohet_JAX_counts,9)

# NIH

NIH_counts_prev <- prevalence_filter(NIH_counts,7)

# Full Dataset
baseline_counts_prev <- prevalence_filter(baseline_counts,20)

## Calculate RS Jensen Shannon distance matrix -- 


# JAX 
JAX.dist <- calculate_rsjensen(JAX_counts_prev)
set.seed(11)
data.adonis=adonis(JAX.dist ~ Sex + Genotype, data=JAX_meta, permutations=10000)
data.adonis$aov.tab 

nohet_JAX.dist <- calculate_rsjensen(nohet_JAX_counts_prev)

data.dist<-nohet_JAX.dist
metadata <- nohet_JAX_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(nohet_JAX.dist ~ Sex + Genotype, data=nohet_JAX_meta, permutations=10000)
data.adonis$aov.tab 


# NIH
NIH.dist <- calculate_rsjensen(NIH_counts_prev)

data.adonis=adonis(NIH.dist ~ Sex + Genotype, data=NIH_meta, permutations=10000)
data.adonis$aov.tab 

# Baseline
baseline.dist <- calculate_rsjensen(baseline_counts_prev)

target <- row.names(baseline.dist)
baseline_meta= baseline_meta[match(target, row.names(baseline_meta)),]
target == row.names(baseline_meta)
baseline.dist <- as.dist(as(baseline.dist, "matrix"))

data.adonis=adonis(baseline.dist ~ Background + Sex + Genotype, data=baseline_meta, permutations=10000)
data.adonis$aov.tab 

# data.adonis=adonis(baseline.dist ~Genotype, data=baseline_meta, permutations=10000)
data.adonis$aov.tab 

## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

nohet_jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=nohet_JAX.dist,
                                         counts = nohet_JAX_counts_prev,
                                         metadata = nohet_JAX_meta,
                                         title="JAX - Fecal Pellet RS Jensen",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = "Baseline/beta_diversity/nohet_JAX_Top_Taxa_PcoA.csv")


jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=JAX.dist,
                                   counts = JAX_counts_prev,
                                   metadata = JAX_meta,
                                   title="JAX - Fecal Pellet RS Jensen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv")

NIH_baseline_pcoa <- generate_pcoA_plots(distance_matrix=NIH.dist,
                                         counts = NIH_counts_prev,
                                         metadata = NIH_meta,
                                         title="NIH - Fecal Pellet RS Jensen",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = "Baseline/beta_diversity/NIH_Top_Taxa_PcoA.csv")

full_baseline_pcoa <- generate_pcoA_plots(distance_matrix=baseline.dist,
                                   counts = baseline_counts_prev,
                                   metadata = baseline_meta,
                                   title="Full Baseline - RS Jensen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Baseline/beta_diversity/Baseline_Top_Taxa_PcoA.csv")

cols <- c("MPTP"="red", "PFF"="black")
generate_pcoA_plots(distance_matrix=JAX.dist,
                    counts = JAX_counts,
                    metadata = JAX_meta,
                    title="JAX - Fecal Pellet RS Jensen",
                    colorvariable = Line,
                    colorvector = cols,
                    wa_scores_filepath = "Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv")

cols <- c("Male"="red", "Female"="black")
generate_pcoA_plots(distance_matrix=JAX.dist,
                    counts = JAX_counts,
                    metadata = JAX_meta,
                    title="JAX - Fecal Pellet RS Jensen",
                    colorvariable = Sex,
                    colorvector = cols,
                    wa_scores_filepath = "Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv")

plot_grid(jax_baseline_pcoa,full_baseline_pcoa, NIH_baseline_pcoa,
          labels=c("E","F", "G"),label_size = 20, nrow=1)
