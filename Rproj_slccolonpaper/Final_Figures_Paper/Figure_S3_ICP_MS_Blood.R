library(here)
library(ggplot2)
library(rlang)
library(cowplot)
library(ggpubr)
library(dplyr)
library(tidyr)
#install.packages("gridGraphics")

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_S3_ICP_MS_Blood.R")

### Data Preprocessing ---
df<- readr::read_csv(here("ICPMS/ICP_MS_Blood.csv"))
row.names(df) <- df$SampleID
df <- df %>% select(-c("SampleID"))

# replace all n/a and declare all element columns as numerical
df[df=="n/a"]<-"0"
vector <- names(df)
elements <- vector[1:7]
df <- df %>% mutate_at(c(elements), as.numeric)
str(df)
df$Genotype_Sex <- paste0(df$Genotype,"_",df$Sex)

# Generate a version of a df without outlier samples 
# create detect outlier function

detect_outlier <- function(x) {
  
  # calculate first quantile
  Quantile1 <- quantile(x, probs=.25,na.rm = TRUE)
  # calculate third quantile
  Quantile3 <- quantile(x, probs=.75,na.rm=TRUE)
  # calculate inter quartile range
  IQR = Quantile3-Quantile1
  # return true or false
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}

# create remove outlier function
remove_outlier <- function(dataframe,
                           columns=names(dataframe)) {
  
  # for loop to traverse in columns vector
  for (col in columns) {
    
    # remove observation if it satisfies outlier function
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  
  # return dataframe
  print("Remove outliers")
  print(dataframe)
  
}


### Figures ---
generate_violin_plots <- function (input_data, X, titlestring) {
  # read in file
  data<-as.data.frame(input_data)
  
  #Ensure correct ordering of levels 
  data$Genotype <- factor(data$Genotype, levels = c("WT","MUT"))
  data$Genotype_Sex <- plyr::revalue(data$Genotype_Sex,c("WT_Male"="WT_M" ,
                                                         "WT_Female"="WT_F",
                                                         "MUT_Female"="MUT_F",
                                                         "MUT_Male"="MUT_M"))
  data$Genotype_Sex <- factor(data$Genotype_Sex, levels=c("WT_M","MUT_M", "WT_F", "MUT_F"))
  
  ggplot(data=data,aes(x={{X}},y=Concentration, fill=Genotype)) + 
    geom_boxplot(alpha=0.25,position=position_dodge(width=.75),size=1,color="black")+
    #scale_shape_manual(values=c(16,10))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25, height=0),alpha=0.8)+
    theme_cowplot(16) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(titlestring)+
    ylab("ug/g") +
    xlab("")
  #ylim(min,max) +
  
}


blood_plots <- list()


elementvector <- c("Iron", "Cobalt", "Copper", "Zinc", "Cadmium", "Manganese", "Selenium")

for (int in 1:7){
  print(int)
  element <- elementvector[int]
  blood <- df %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)
 
  blood_plots[[int]] <- blood

}

## Assemble Multi-Panel Figure
cd_blood <- blood_plots[[5]] + ylim(0,0.05)

plot_grid(blood_plots[[1]],blood_plots[[2]],
                     blood_plots[[3]],blood_plots[[4]],
                     cd_blood,blood_plots[[6]],blood_plots[[7]],
                      nrow = 2, ncol=4,
                      labels=c("A", "","","","",""
                      ),label_size = 20) 

fig_2_bottom <- plot_grid(blood_plots[[1]],blood_plots[[2]],
          blood_plots[[3]],blood_plots[[4]],
          cd_blood,blood_plots[[6]],
          nrow = 1, ncol=6,
          labels=c("D", "","","","",""
          ),label_size = 20) 
dev.new(width=10,height=10)
fig_2_bottom
### Figures by Sex_Genotype---


blood_plots <- list()


elementvector <- c("Iron", "Cobalt", "Copper", "Zinc", "Cadmium", "Manganese", "Selenium")

for (int in 1:7){
  print(int)
  element <- elementvector[int]
  blood <- df %>% 
    remove_outlier(element) %>%
    pivot_longer(cols=all_of(elementvector),
                 names_to="Element",
                 values_to="Concentration") %>%
    filter(Element==element) %>%
    generate_violin_plots(X=Genotype, titlestring=element)+ facet_wrap(~Sex)
  
  blood_plots[[int]] <- blood
  
}

## Assemble Multi-Panel Figure
cd_blood <- blood_plots[[5]] + ylim(0,0.05)

plot_grid(blood_plots[[1]],blood_plots[[2]],
          blood_plots[[3]],blood_plots[[4]],
          cd_blood,blood_plots[[6]],blood_plots[[7]],
          nrow = 1, ncol=6,
          labels=c("A", "","","","",""
          ),label_size = 20) 

fig_s3_bottom <- plot_grid(blood_plots[[1]],blood_plots[[2]],
                        blood_plots[[3]],blood_plots[[4]],
                        cd_blood,blood_plots[[6]],
                        nrow = 1, ncol=6,
                        labels=c("D", "","","","",""
                        ),label_size = 20) 

dev.new(width=10,height=10)
fig_s3_bottom

### Statistics On Full Dataset--- 

df<- readr::read_csv(here("ICPMS/ICP_MS_Blood.csv"))
row.names(df) <- df$SampleID
df <- df %>% select(-c("SampleID"))

# replace all n/a and declare all element columns as numerical
df[df=="n/a"]<-"0"
vector <- names(df)
elements <- vector[1:7]
df <- df %>% mutate_at(c(elements), as.numeric)
str(df)
df$Genotype_Sex <- paste0(df$Genotype,"_",df$Sex)

element_stats_para <- list()
element_stats_nonpara <-list()

#filter to just males or just females
df_f <- df %>% filter(Sex=="Female")
df_m <- df %>% filter(Sex=="Male")


for (int in 1:7){
  print(int)
  bld_para <- t.test(as.numeric(unlist(df_f[,int]))~Genotype,df_f)
  bld_nonpara <- wilcox.test(as.numeric(unlist(df_f[,int]))~Genotype,df_f)

  element_stats_para[[int]] <-list(print(bld_para))
  element_stats_nonpara[[int]] <-list(print(bld_nonpara))
  
}


# Iron 
element_stats_para[[1]]
element_stats_nonpara[[1]]

# Cobalt 
element_stats_para[[2]]
element_stats_nonpara[[2]]

# Copper
element_stats_para[[3]]
element_stats_nonpara[[3]]

# Zinc
element_stats_para[[4]]
element_stats_nonpara[[4]]

# Cadmium
element_stats_para[[5]]
element_stats_nonpara[[5]]

# Manganese
element_stats_para[[6]]
element_stats_nonpara[[6]]

# Selenium
element_stats_para[[7]]
element_stats_nonpara[[7]]

for (int in 1:7){
  print(int)
  bld_para <- t.test(as.numeric(unlist(df_m[,int]))~Genotype,df_m)
  bld_nonpara <- wilcox.test(as.numeric(unlist(df_m[,int]))~Genotype,df_m)
  
  element_stats_para[[int]] <-list(print(bld_para))
  element_stats_nonpara[[int]] <-list(print(bld_nonpara))
  
}


# Iron 
element_stats_para[[1]]
element_stats_nonpara[[1]]

# Cobalt 
element_stats_para[[2]]
element_stats_nonpara[[2]]

# Copper
element_stats_para[[3]]
element_stats_nonpara[[3]]

# Zinc
element_stats_para[[4]]
element_stats_nonpara[[4]]

# Cadmium
element_stats_para[[5]]
element_stats_nonpara[[5]]

# Manganese
element_stats_para[[6]]
element_stats_nonpara[[6]]

# Selenium
element_stats_para[[7]]
element_stats_nonpara[[7]]

