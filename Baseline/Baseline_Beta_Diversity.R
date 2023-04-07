library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/SLC_GitHub/slcproject/SLC_Microbiome_Baseline/beta_diversity/")

generate_cs_pcoA_plots <- function(ordination_file, metadata, title, colorvariable,colorvector){
  data<-read.csv(ordination_file, header=FALSE)
 # data<-read.csv("rpca/pcoa_rpca_s5_min10000_PFF_Jejunum_min10000_no_tax_PFF_ASV_table.qza.txt/ordination.csv", header=FALSE)
  
  metadata <- read.delim(metadata, header=TRUE)
  #metadata <- read.delim("../starting_files/PFF_Mapping.tsv", header=TRUE)
  #metadata$SampleID <- metadata$X.SampleID
  
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
  data$Genotype<-factor(data$Genotype, levels=c("WT", "HET", "MUT"))
  #data <- data %>% arrange(Timepoint) 
  
  p<- ggplot(data, aes(x=PC1, y=PC2, colour={{colorvariable}})) + 
    geom_point(size=3) + 
    scale_colour_manual(name="",values={{colorvector}}) +
    #scale_shape_manual(name="", values=c(16,10)) +
    #scale_color_viridis_d()+
    xlab(str_PC1) +
    ylab(str_PC2) +
    theme_cowplot(16)+
    #stat_ellipse()+
    theme(legend.position="top",legend.justification = "center") +
    #geom_line(aes(group = MouseID),color="darkgrey", arrow = arrow(type = "closed",length=unit(0.075, "inches")))+
    #geom_point(aes(x = PC1, y = PC2, shape = Timepoint), size = 3) + 
    #geom_path(aes(x = PC1, y = PC2, group = MouseID), arrow = arrow(length = unit(0.55, "cm")))
  #coord_fixed(ratio=1/2)+
  labs(title= paste0({{title}})) 
  p
}

genotype_cols<- c("MUT" = "red", "WT" = "black", "HET"="blue")

baseline_rpca<- generate_cs_pcoA_plots("pcoa_rpca_s20_min10000_Baseline_ASV_table_Silva_v138_1.qza.txt/ordination.csv", "../starting_files/Baseline_Metadata - Baseline_Metadata.tsv", "Baseline RPCA", Genotype,genotype_cols)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse()
baseline_rpca

baseline_rpca_by_Line <- generate_cs_pcoA_plots("pcoa_rpca_s20_min10000_Baseline_ASV_table_Silva_v138_1.qza.txt/ordination.csv", "../starting_files/Baseline_Metadata - Baseline_Metadata.tsv", "Baseline RPCA", Genotype,genotype_cols)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse()+
  facet_wrap(~Line)
baseline_rpca_by_Line

baseline_rpca_by_Sex <- generate_cs_pcoA_plots("pcoa_rpca_s20_min10000_Baseline_ASV_table_Silva_v138_1.qza.txt/ordination.csv", "../starting_files/Baseline_Metadata - Baseline_Metadata.tsv", "Baseline RPCA", Genotype,genotype_cols)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse()+
  facet_wrap(~Sex)
baseline_rpca_by_Sex

plot_grid(baseline_rpca, baseline_rpca_by_Line, baseline_rpca_by_Sex, labels=c("A","B", "C"))

baseline_bc<- generate_cs_pcoA_plots("bray_curtis/s20_min10000_Baseline_ASV_table_Silva_v138_1_pcoa/ordination.csv", "../starting_files/Baseline_Metadata - Baseline_Metadata.tsv", "Baseline bc", Genotype,genotype_cols)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse()
baseline_bc

baseline_bc_by_Line<- generate_cs_pcoA_plots("bray_curtis/s20_min10000_Baseline_ASV_table_Silva_v138_1_pcoa/ordination.csv", "../starting_files/Baseline_Metadata - Baseline_Metadata.tsv", "Baseline bc", Genotype,genotype_cols)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse()+
  facet_wrap(~Line)
baseline_bc_by_Line

baseline_bc_by_Sex<- generate_cs_pcoA_plots("bray_curtis/s20_min10000_Baseline_ASV_table_Silva_v138_1_pcoa/ordination.csv", "../starting_files/Baseline_Metadata - Baseline_Metadata.tsv", "Baseline bc", Genotype,genotype_cols)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse()+
  facet_wrap(~Sex)
baseline_bc_by_Sex

plot_grid(baseline_bc, baseline_bc_by_Line, baseline_bc_by_Sex, labels=c("A","B", "C"))

### Beta - Diversity Stats ---

## RPCA --
metadata <-read.delim("../starting_files/Baseline_Metadata - Baseline_Metadata.tsv",sep="\t",header=TRUE, row.names=1) #mapping file
write.csv(metadata, "../starting_files/Baseline_Metadata.csv")

data.dist<-read.table(file= "dm_rpca_s20_min10000_Baseline_ASV_table_Silva_v138_1.qza.txt/distance-matrix.tsv")
metadata <- read.csv("../starting_files/Baseline_Metadata.csv", header=TRUE, row.names=1)

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~Line+ Sex+ Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

data.adonis=adonis(data.dist ~ Line + Sex*Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

## RPCA --
metadata <-read.delim("../starting_files/Baseline_Metadata - Baseline_Metadata.tsv",sep="\t",header=TRUE, row.names=1) #mapping file
write.csv(metadata, "../starting_files/Baseline_Metadata.csv")

data.dist<-read.table(file= "bray_curtis/export_s20_min10000_Baseline_ASV_table_Silva_v138_1_DM/distance-matrix.tsv")
metadata <- read.csv("../starting_files/Baseline_Metadata.csv", header=TRUE, row.names=1)

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~Line+ Sex+ Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

data.adonis=adonis(data.dist ~ Line + Sex*Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

data.adonis=adonis(data.dist ~ Sex*Line + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab
