library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Long_Term_RS_Jensen_Shannon.R")

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
lumcol_counts<- lumcol_counts%>% filter(prevalence>=13) #[Result: 557 features, 90 samples]

lumcol_counts <- lumcol_counts %>% select(-c(prevalence))

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
muccol_counts <- muccol_counts %>% filter(prevalence>=12) #[Result: 596 features, 89 samples]

muccol_counts <- muccol_counts %>% select(-c(prevalence))

## Calculate RS Jensen Shannon distance matrix -- 


muccol.dist <- calculate_rsjensen(muccol_counts)
lumcol.dist <- calculate_rsjensen(lumcol_counts)


## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lumcol.dist,
                                     counts = lumcol_counts,
                                     metadata = lumcol_meta,
                                     title="Long Term - Luminal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Long_Term/LumCol_Top_Taxa_PcoA.csv")

slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                     counts = muccol_counts,
                                     metadata = muccol_meta,
                                     title="Long Term - Mucosal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Long_Term/MucCol_Top_Taxa_PcoA.csv")

plot_grid(slt_lc_pcoa,slt_mc_pcoa, labels=c("C","D"),label_size = 20)
