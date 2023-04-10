```R
> # Iron 
> element_stats_para[[1]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -2.7287, df = 37.376, p-value = 0.009636
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -234.8644  -34.7384
sample estimates:
mean in group MUT  mean in group WT 
         677.1142          811.9157 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = -0.41903, df = 41.212, p-value = 0.6774
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -138.48442   90.88552
sample estimates:
mean in group MUT  mean in group WT 
         366.5477          390.3472 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = 1.1198, df = 32.461, p-value = 0.271
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -76.67843 264.17162
sample estimates:
mean in group MUT  mean in group WT 
         328.1754          234.4288 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = 0.46633, df = 43.904, p-value = 0.6433
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -45.06846  72.20123
sample estimates:
mean in group MUT  mean in group WT 
         189.8863          176.3200 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 1.2857, df = 20.054, p-value = 0.2132
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -32.5599 137.2244
sample estimates:
mean in group MUT  mean in group WT 
         180.7050          128.3727 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = 1.2093, df = 37.605, p-value = 0.2341
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -8.790092 34.849392
sample estimates:
mean in group MUT  mean in group WT 
         98.25350          85.22385 


> element_stats_nonpara[[1]]
[[1]]

	Wilcoxon rank sum exact test

data:  df_fp_col[, int] by Genotype
W = 113, p-value = 0.01074
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum exact test

data:  df_fp_si[, int] by Genotype
W = 235, p-value = 0.6939
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_col[, int] by Genotype
W = 262, p-value = 0.2771
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum exact test

data:  df_muc_si[, int] by Genotype
W = 307, p-value = 0.6976
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum exact test

data:  df_ts_col[, int] by Genotype
W = 267, p-value = 0.07178
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum exact test

data:  df_ts_si[, int] by Genotype
W = 242, p-value = 0.2648
alternative hypothesis: true location shift is not equal to 0
```