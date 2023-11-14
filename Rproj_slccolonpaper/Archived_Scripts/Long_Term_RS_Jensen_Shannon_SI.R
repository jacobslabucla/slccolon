library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Long_Term_RS_Jensen_Shannon.R")

metadata <- read.table("Long_Term/starting_files/SLC_LT_metadata.tsv", header=TRUE)
counts <- read.table("Long_Term/starting_files/SLT_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal SI
lumsi_meta <- metadata %>% filter(Subset=="Luminal_SI", SampleID %in% names(counts))
row.names(lumsi_meta) <- lumsi_meta$SampleID
lumsi <- lumsi_meta$SampleID
lumsi_counts <- counts %>% select(all_of(lumsi))

# Mucosal SI
mucsi_meta <- metadata %>% filter(Subset=="Mucosal_SI", SampleID %in% names(counts))
row.names(mucsi_meta) <- mucsi_meta$SampleID
mucsi <- mucsi_meta$SampleID
mucsi_counts <- counts %>% select(all_of(mucsi))

# Luminal Colon no HET 
nohet_lt_lumsi_meta <- metadata %>%
  filter(Genotype!="HET")%>%
  filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(nohet_lt_lumsi_meta) <- nohet_lt_lumsi_meta$SampleID
lumsi <- nohet_lt_lumsi_meta$SampleID
nohet_lumsi_counts <- counts %>% select(all_of(lumsi))

# Mucosal Colon no HET
nohet_lt_mucsi_meta <- metadata %>% 
  filter(Genotype!="HET")%>%
  filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(nohet_lt_mucsi_meta) <- nohet_lt_mucsi_meta$SampleID
mucsi <- nohet_lt_mucsi_meta$SampleID
nohet_mucsi_counts <- counts %>% select(all_of(mucsi))

## Prevalence filter datasets -- 
# Luminal SI
0.15*90 #13 samples
lumsi_counts <- prevalence_filter(lumsi_counts,13)

# Mucosal SI
0.15*89
mucsi_counts <- prevalence_filter(mucsi_counts,13)

# Luminal Colon - no HET 
0.15*42 #6 samples
nohet_lt_lumsi_counts <- prevalence_filter(nohet_lumsi_counts,6)

# Mucosal Colon - no HET 
0.15*41 #6 samples
nohet_lt_mucsi_counts <- prevalence_filter(nohet_mucsi_counts,6)

## Calculate RS Jensen Shannon distance matrix -- 


mucsi.dist <- calculate_rsjensen(mucsi_counts)
lumsi.dist <- calculate_rsjensen(lumsi_counts)
lt_lumsi.dist <- calculate_rsjensen(nohet_lt_lumsi_counts)
lt_mucsi.dist <- calculate_rsjensen(nohet_lt_mucsi_counts)

## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

slt_lsi_pcoa <- generate_pcoA_plots(distance_matrix=lumsi.dist,
                                   counts = lumsi_counts,
                                   metadata = lumsi_meta,
                                   title=" Luminal SI",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/lumsi_Top_Taxa_PcoA.csv")

slt_msi_pcoa <- generate_pcoA_plots(distance_matrix=mucsi.dist,
                                   counts = mucsi_counts,
                                   metadata = mucsi_meta,
                                   title="Mucosal SI",
                                   colorvariable = Genotype,
                                   colorvector = cols,
                                   wa_scores_filepath = "Long_Term/mucsi_Top_Taxa_PcoA.csv")
nohet_slt_lc_pcoa <- generate_pcoA_plots(distance_matrix=lt_lumsi.dist,
                                         counts = nohet_lumsi_counts,
                                         metadata = nohet_lt_lumsi_meta,
                                         title="Long Term - Luminal Colon RS Jensen",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = "Long_Term/nohet_lumsi_Top_Taxa_PcoA.csv")

nohet_slt_mc_pcoa <- generate_pcoA_plots(distance_matrix=lt_mucsi.dist,
                                         counts = nohet_mucsi_counts,
                                         metadata = nohet_lt_mucsi_meta,
                                         title="Long Term - Mucosal Colon RS Jensen",
                                         colorvariable = Genotype,
                                         colorvector = cols,
                                         wa_scores_filepath = "Long_Term/nohet_mucsi_Top_Taxa_PcoA.csv")

plot_grid(slt_lsi_pcoa,slt_msi_pcoa, labels=c("C","D"),label_size = 20)
plot_grid(nohet_slt_lc_pcoa,nohet_slt_mc_pcoa, labels=c("C","D"),label_size = 20)


## PERMANOVA

# Luminal SI
data.dist<-lumsi.dist
metadata <- lumsi_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal SI
data.dist<-mucsi.dist
metadata <- mucsi_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Luminal Colon -- no HET 
data.dist<-lt_lumsi.dist
metadata <- nohet_lt_lumsi_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon -- no HET
data.dist<-lt_mucsi.dist
metadata <- nohet_lt_mucsi_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab


## Repeat-Measures-Aware 
# Luminal Colon
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(lumsi.dist,
                       lumsi_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(lt_lumsi.dist,
                       nohet_lt_lumsi_meta,
                       site,
                       mouseID,
                       order_vector)

# Mucosal Colon
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(mucsi.dist,
                       mucsi_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sex","Genotype","MouseID")
order_vector <- c("Sex","Site","Genotype")
run_repeated_PERMANOVA(lt_mucsi.dist,
                       nohet_lt_mucsi_meta,
                       site,
                       mouseID,
                       order_vector)