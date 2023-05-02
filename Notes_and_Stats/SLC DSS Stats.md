### Histology
##### MUT vs WT 
```R
> df_mut <- histology %>% filter(Genotype!="HET")
> wilcox.test(Score~Genotype,df_mut)

	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 33.5, p-value = 0.9191
alternative hypothesis: true location shift is not equal to 0
```

##### HET vs WT 
```R
> df_het <- histology %>% filter(Genotype!="MUT")
> wilcox.test(Score~Genotype, df_het)

	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 28, p-value = 0.6916
alternative hypothesis: true location shift is not equal to 0
```