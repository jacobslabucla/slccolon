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

# full dataset
baseline_meta <- metadata
row.names(baseline_meta)=baseline_meta$SampleID
baseline_counts <- counts 

## Prevalence filter datasets -- 
# JAX
t_df_input_data<-as.data.frame(t(JAX_counts))

ctr= 0
prevalence <- vector(mode="numeric")

for(i in 1:ncol(t_df_input_data)){
  v<-t_df_input_data %>% pull(i)
  for(j in 1:length(v)){
    if (t_df_input_data[j,i]>0){
      ctr=1+ctr
    }
    else {
      ctr=ctr
    }
  }
  prevalence<-append(prevalence,ctr)
  ctr=0
}
0.15*87 #13.5 samples 
JAX_counts$prevalence<-prevalence #features present in at least 13 samples our of 90
JAX_counts_prev<- JAX_counts%>% filter(prevalence>=13) #[Result: 557 features, 90 samples]

JAX_counts_prev <- JAX_counts %>% select(-c(prevalence))

# Full Dataset
t_df_input_data<-as.data.frame(t(baseline_counts))

ctr= 0
prevalence <- vector(mode="numeric")

for(i in 1:ncol(t_df_input_data)){
  v<-t_df_input_data %>% pull(i)
  for(j in 1:length(v)){
    if (t_df_input_data[j,i]>0){
      ctr=1+ctr
    }
    else {
      ctr=ctr
    }
  }
  prevalence<-append(prevalence,ctr)
  ctr=0
}
0.15*133 #20 samples 
baseline_counts$prevalence<-prevalence #features present in at least 15 samples our of 100
baseline_counts_prev <- baseline_counts %>% filter(prevalence>=20) #[Result: 244 features, 133 samples]

baseline_counts_prev <- baseline_counts %>% select(-c(prevalence))

## Calculate RS Jensen Shannon distance matrix -- 


baseline.dist <- calculate_rsjensen(baseline_counts_prev)
JAX.dist <- calculate_rsjensen(JAX_counts_prev)

data.adonis=adonis(JAX.dist ~ Sex + Genotype, data=JAX_meta, permutations=10000)
data.adonis$aov.tab 

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

jax_baseline_pcoa <- generate_pcoA_plots(distance_matrix=JAX.dist,
                                   counts = JAX_counts,
                                   metadata = JAX_meta,
                                   title="JAX - Fecal Pellet RS Jensen",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Baseline/beta_diversity/JAX_Top_Taxa_PcoA.csv")

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

plot_grid(jax_baseline_pcoa,full_baseline_pcoa, labels=c("E","F"),label_size = 20)
