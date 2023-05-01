library(ggplot2)

here::i_am("Rproj_slccolonpaper/Overlapping_Genera.R")

## Luminal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Luminal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
print(data$feature)

trios_l6_het_significant <- data %>% filter(value=="HET")
trios_l6_mut_significant <- data %>% filter(value=="MUT")

make_genus_level_taxa_dotplot(trios_l6_mut_significant, 
                              "Trios/differential_taxa/L6_Luminal_ColonRelative_Abundance-ASV.RDS",
                              "Trios Luminal Colon MUT: Genus ~ SeqRun + Site + Sex + Genotype + (1|MouseID)")
phylum_colors = c("Red","blue","green")
## Make a taxonomy dotplot --

make_genus_level_taxa_dotplot <- function(ASV_significant_results_dataset,
                              Relative_Abundance_filepath_rds,titlestring){
  data <- data %>% filter(qval <0.25)
  data <- data %>% filter(metadata=="Genotype")
  data$Taxon <- data$feature
  data$Phylum <- gsub(".*p__","",data$Taxon)
  data$Phylum <- gsub("\\..*","",data$Phylum)
  data$Family<- gsub(".*f__","",data$Taxon)
  data$Family <-  gsub("\\..*","",data$Family)
  data$Order<- gsub(".*o__","",data$Taxon)
  data$Order <-  gsub("\\..*","",data$Order)
  data$Genus<- gsub(".*g__","",data$Taxon)
  
  data$annotation <- gsub("\\.E","E",data$Genus)
  data$annotation <- gsub("\\.","_",data$annotation)
  data$annotation <- gsub("__","_",data$annotation)
  #data$Genus <- gsub("\\..*","",data$Genus)
  data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
  data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))
  
  #append relative abundance data 
  #relA <- readRDS("Trios/differential_taxa/L6_Luminal_ColonRelative_Abundance-ASV.RDS")
  relA <- readRDS(Relative_Abundance_filepath_rds)
  relA$feature <- row.names(relA)
  relA$feature <- gsub(";",".",relA$feature)
  relA$feature <- gsub(" ",".",relA$feature)
  relA$Relative_Abundance <- relA$V1
  data<-merge(data,relA,by="feature")
  print(data$feature)
  print(summary(data$Relative_Abundance))
  max(data$Relative_Abundance)
  
  #make graph
  y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  data$annotation= factor(as.character(data$annotation), levels = names(y))
  baseline_DAT <- ggplot(data, aes(x = coef, y = annotation ,color=Phylum)) + 
    geom_point(aes(size = sqrt(Relative_Abundance))) + 
    scale_size_continuous(name="Relative Abundance",range = c(0.5,8),
                          limits=c(sqrt(0.000001),sqrt(0.3)),
                          breaks=c(sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),
                          labels=c("0.0001","0.001","0.01","0.1")) + 
    #scale_color_manual(name="Phylum", values = phyla_colors)+
    geom_vline(xintercept = 0) + 
    xlab(label="Log2 Fold Change")+
    ylab(label=NULL)+
    theme_cowplot(16) +
    ggtitle(titlestring) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "right") 
  baseline_DAT 
  
}
