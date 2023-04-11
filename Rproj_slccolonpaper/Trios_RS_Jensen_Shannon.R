library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_ASV_level_Maaslin2.R")

metadata <- read.table("Trios/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

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
lumcol_counts<- lumcol_counts%>% filter(prevalence>=13) #[Result: 383 features, 90 samples]

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
muccol_counts <- muccol_counts %>% filter(prevalence>=12) #[Result: 450 features, 83 samples]

muccol_counts <- muccol_counts %>% select(-c(prevalence))

## Calculate RS Jensen Shannon distance matrix -- 


muccol.dist <- calculate_rsjensen(muccol_counts)
lumcol.dist <- calculate_rsjensen(lumcol_counts)

data.adonis=adonis(lumcol.dist ~ Litter + Site +Sex + Genotype, data=lumcol_meta, permutations=10000)
data.adonis$aov.tab

## Principal Coordinates Analysis -- 
# calculate mds
mds <- cmdscale(lumcol.dist, eig = TRUE, x.ret = TRUE)

mds_values <- mds$points
wa_scores <- wascores(mds_values, t(lumcol_counts))
wa_scores <- data.frame(sample = rownames(wa_scores),
                        x = wa_scores[,1],
                        y = wa_scores[,2])

# isolate taxa with strongest contribution to principal coordinate axes
n_taxa <- 10
wa_scores_1<- head(arrange(wa_scores, desc(abs(wa_scores$x))), n = n_taxa)
wa_scores_2<- head(arrange(wa_scores, desc(abs(wa_scores$y))), n = n_taxa)
wa_scores_final <- rbind(wa_scores_1, wa_scores_2)

# calculate percentage of variation that each mds axis accounts for
mds_var_per <- round(mds$eig/sum(mds$eig) * 100, 1)

# plot
mds_data <- data.frame(sample = rownames(mds_values),
                       x = mds_values[,1],
                       y = mds_values[,2])

#merge phenotypic data 
lumcol_meta$sample <- lumcol_meta$SampleID
mds_meta <- merge(mds_data, lumcol_meta, by = "sample")

mds_plot <- ggplot(mds_meta, aes(x, y, color = Genotype)) +
  geom_point(size = 2, alpha = 0.85) +
 # scale_color_manual(values = za_pal) +
  labs(x = paste("PC1(", mds_var_per[1], "%)",sep=""),
       y = paste("PC2(", mds_var_per[2], "%)",sep=""),
       color = "") +
  theme_cowplot(16) +
  coord_fixed() +
  theme(legend.position = "top",
        legend.justification = "center",
        plot.margin = margin(t = 0, unit = "cm"),
        legend.margin = margin(b = -0.25, unit = "cm"))


generate_pcoA_plots <- function(distance_matrix, metadata, title, colorvariable,colorvector){
  data<-read.csv(ordination_file, header=FALSE)
  metadata <- read.delim(metadata, header=TRUE,row.names = 1)
  metadata$SampleID<-row.names(metadata)
  #store PC1 and Pc2
  PC1<-data[5,1]
  PC1 <-round(as.numeric(PC1)*100, digits=1)
  PC2<-data[5,2]
  PC2 <-round(as.numeric(PC2)*100, digits=1)
  PC1 <-as.character(PC1)
  str_PC1<-paste0("PC 1 (", PC1,"%)")
  str_PC2<-paste0("PC 2 (", PC2, "%)")
  
  #edit dataframe
  data<-data[,1:4]
  data <- slice(data, 1:(n() - 4))     # Apply slice & n functions
  data<-data[-c(1,2,3,4,5,6,7,8,9),]
  
  #rename columns
  names(data)[names(data) == "V1"] <- "SampleID" 
  names(data)[names(data)=="V2"] <- "PC1" 
  names(data)[names(data)=="V3"] <- "PC2"
  names(data)[names(data)=="V4"] <- "PC3"
  # data$SampleID<-gsub(".","",data$SampleID)
  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #declare factors
  data$Diet<-data$Group
  data$Diet<-factor(data$Diet, levels=c("Control", "BD"))
  data$Diet <- plyr::revalue(data$Diet, c("Control"="CD","BD"="BD"))
  
  p<- ggplot(data, aes(x=PC1, y=PC2, colour={{colorvariable}})) + 
    geom_point(size=3) + 
    scale_colour_manual(name="",values={{colorvector}}) +
    #scale_color_viridis_d()+
    xlab(str_PC1) +
    ylab(str_PC2) +
    theme_cowplot(16)+
    theme(legend.position="top",legend.justification = "center") +
    #geom_text(nudge_y = .05) +
    #geom_line(aes(group = MouseID),color="darkgrey", arrow = arrow(type = "closed",length=unit(0.075, "inches")))+
    #geom_point(aes(x = PC1, y = PC2, shape = Timepoint), size = 3) + 
    #geom_path(aes(x = PC1, y = PC2, group = MouseID), arrow = arrow(length = unit(0.55, "cm")))+
    #coord_fixed(ratio=1/2)+
    labs(title= paste0({{title}})) 
  p
}
