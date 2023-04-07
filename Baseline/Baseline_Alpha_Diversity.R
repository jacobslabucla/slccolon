
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(ggpubr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/SLC_GitHub/slcproject/SLC_Microbiome_Baseline/")

otus<- readr::read_delim("alpha_diversity/alpha_min_10000_Baseline_ASV_table_Silva_v138_1/otus_dir/alpha-diversity.tsv")
row.names(otus) <- otus$...1
shannon<-readr::read_delim("alpha_diversity/alpha_min_10000_Baseline_ASV_table_Silva_v138_1/shannon_dir/alpha-diversity.tsv")
row.names(shannon) <- shannon$...1

data<- merge(otus,shannon, by="...1")
data$SampleID <- data$...1

metadata<- read.delim("starting_files/Baseline_Metadata - Baseline_Metadata.tsv", row.names=1)
metadata$SampleID <- row.names(metadata)
data_meta <- merge(data,metadata, by="SampleID")

### Split into stool and intestine datasets
### Function for plotting alpha diversity ---
generate_adiv_plots <- function(input_data, X, Y, min, max){
  #read in files
  data <- as.data.frame(input_data)

  #declare order of variables
  data$Genotype <- factor(data$Genotype, levels=c("WT", "HET","MUT"))
  #graph plot
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill=Genotype)) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    #geom_line(aes(group = MouseID,color=Genotype),size=1)+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(16) +
    ylim(min,max) +
    theme(legend.position = "none")
  
}

### Make and store plots ---
compare <-c(c("WT","HET"), c("WT","MUT"))

adiv_baseline_shannon<- generate_adiv_plots(data_meta, Genotype, shannon, 2, 7) +
  facet_grid(~Line)+
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("Shannon")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_baseline_shannon

adiv_baseline_otus<- generate_adiv_plots(data_meta, Genotype, observed_otus, 0, 300) +
  facet_grid(~Line)+
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("# ASVs")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_baseline_otus

plot_grid(adiv_baseline_shannon, adiv_baseline_otus, labels=c("A","B"))
### Alpha Diversity Stats ---
#stool
stool$Genotype <-factor(stool$Genotype, levels=c("WT", "KO"))
output <- lme(fixed= shannon ~ Timepoint*Genotype, random = ~1|MouseID, data=stool)
summary(output)
output <- lme(fixed= observed_otus ~ Timepoint*Genotype, random = ~1|MouseID, data=stool)
summary(output)
output <- lme(fixed= shannon ~ Timepoint+Genotype, random = ~1|MouseID, data=stool)
summary(output)
output <- lme(fixed= observed_otus ~ Timepoint+Genotype, random = ~1|MouseID, data=stool)
summary(output)

stool_day0 <- stool %>% filter(Timepoint=="Day0")
stool_day7 <- stool %>% filter(Timepoint=="Day7 ")

wilcox.test(shannon~Genotype,stool_day0)
wilcox.test(shannon~Genotype,stool_day7)
wilcox.test(observed_otus~Genotype,stool_day0)
wilcox.test(observed_otus~Genotype,stool_day7)

t.test(shannon~Genotype,stool_day0)
t.test(shannon~Genotype,stool_day7)
t.test(observed_otus~Genotype,stool_day0)
t.test(observed_otus~Genotype,stool_day7)


#intestine
t.test(shannon~Genotype,intestine)
t.test(observed_otus~Genotype,intestine)
wilcox.test(observed_otus~Genotype,intestine)
wilcox.test(shannon~Genotype,intestine)
