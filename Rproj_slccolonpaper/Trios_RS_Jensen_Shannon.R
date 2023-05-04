library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_RS_Jensen_Shannon.R")

metadata <- read.table("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

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
t_df_input_data<-as.data.frame(t(lumcol_counts))

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
0.15*90 #13.5 samples 
lumcol_counts$prevalence<-prevalence #features present in at least 13 samples our of 90
lumcol_counts<- lumcol_counts%>% filter(prevalence>=13) #[Result: 383 features, 90 samples]

trios_lumcol_counts <- lumcol_counts %>% select(-c(prevalence))

# Mucosal Colon 
t_df_input_data<-as.data.frame(t(muccol_counts))

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
0.15*82 #12 samples 
muccol_counts$prevalence<-prevalence #features present in at least 15 samples our of 100
muccol_counts <- muccol_counts %>% filter(prevalence>=12) #[Result: 450 features, 83 samples]

trios_muccol_counts <- muccol_counts %>% select(-c(prevalence))

## Calculate RS Jensen Shannon distance matrix -- 


trios_muccol.dist <- calculate_rsjensen(trios_muccol_counts)
trios_lumcol.dist <- calculate_rsjensen(trios_lumcol_counts)



## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=lumcol.dist,
                    counts = trios_lumcol_counts,
                    metadata = trios_lumcol_meta,
                    title="Trios - Luminal Colon RS Jensen",
                    colorvariable = Genotype,
                    colorvector = cols,
                    wa_scores_filepath = "Trios/beta_diversity/LumCol_Top_Taxa_PcoA.csv")

trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                     counts = trios_muccol_counts,
                                     metadata = trios_muccol_meta,
                                     title="Trios - Mucosal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/MucCol_Top_Taxa_PcoA.csv")

plot_grid(trios_lc_pcoa, trios_mc_pcoa,labels=c("A","B"),label_size=20)
