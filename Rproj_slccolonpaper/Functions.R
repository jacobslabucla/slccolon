## Make a taxa summary plot --
generate_L6_taxa_plots <- function(path_to_RDS, titlestring,greppattern, fillvector){
  #L2_lum<-readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/SI_LumMuc_L6.RDS")
  #taxa <- gsub(".*g__","",taxa)
  #cols<-assign_cols
  titlestring<-c(titlestring)
  L2_lum<-readRDS(path_to_RDS)
  L2_lum<- as.matrix(L2_lum)
  L2_lum<-funrar::make_relative(L2_lum)
  L2_lum<-as.data.frame(t(L2_lum))
  toptaxa<- rowMeans(L2_lum)
  L2_lum$averageRA <-toptaxa
  L2_lum <- L2_lum %>% dplyr::mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
  L2_lum <-select(L2_lum,-averageRA)
  
  taxa<-L2_lum$keeptaxa
  L2_lum <- select(L2_lum,-keeptaxa)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  taxa<-gsub(greppattern,"",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  L2_lum$Site <- factor(L2_lum$Site, levels=c("WT", "HET","MUT"))
  cols <- fillvector
  ggplot2::ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    theme(legend.position = "none")+
    theme_cowplot(20) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

## handle the genera names for taxonomy boxplot --
wrangle_genera_names <- function(csv_dataframe, filepathstring, rds_string){
  input_data<-read.csv(csv_dataframe, row.names=1, header=TRUE)
  taxa<-colnames(input_data)
  colnames <- strsplit(taxa, ".f__")
  
  family=rlang::new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(colnames)) {
    family[i] <- colnames[[i]][2]
    i=i+1
  }
  
  family<-unlist(family)
  family <- strsplit(family, ".g__")
  
  genus =new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(family)) {
    genus[i] <- family[[i]][2]
    i=i+1
  }
  
  family<-as.list(family)
  
  i=1
  for (i in 1:length(genus)) {
    if (isFALSE(genus[[i]]=="NA")) {
      genus[[i]] = genus[[i]] 
    }
    else {
      
      genus[[i]] <- paste0(family[[i]]," (f)")   
    }
    i=i+1
  }
  colnames(input_data) <-as.character(genus)
  
  saveRDS(input_data, paste0(filepathstring, rds_string))
  
}

## Make a taxonomy dotplot --

make_taxa_dotplot <- function(ASV_significant_results_dataset, taxonomy_tsv_filepath, 
                              Relative_Abundance_filepath_rds,titlestring){
data<-as.data.frame(ASV_significant_results_dataset)
taxonomy <- read.delim(taxonomy_tsv_filepath)
taxonomy$feature <- taxonomy$Feature.ID
data <- merge(data,taxonomy, by="feature")
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub(";.*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub(";.*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub(";.*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)
data$Genus <-  gsub(";.*","",data$Genus)
data$Species <- gsub(".*s__","",data$Taxon)
data$annotation <- paste0(data$Genus," ", data$Species)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

#append relative abundance data 
relA <- readRDS(Relative_Abundance_filepath_rds)
relA$feature <- row.names(relA)
relA$Relative_Abundance <- relA$V1
data<-merge(data,relA,by="feature")
min(data$Relative_Abundance)
max(data$Relative_Abundance)

#make graph
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
baseline_DAT <- ggplot(data, aes(x = coef, y = annotation, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Relative Abundance",range = c(0.5,8),
                        limits=c(sqrt(0.00001),sqrt(0.3)),
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
