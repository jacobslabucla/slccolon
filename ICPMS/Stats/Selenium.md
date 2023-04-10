```R
> # Selenium
> element_stats_para[[7]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -1.2878, df = 37.366, p-value = 0.2057
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4681808  0.1042475
sample estimates:
mean in group MUT  mean in group WT 
         1.670700          1.852667 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = -0.40987, df = 37.082, p-value = 0.6843
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.5862924  0.3889920
sample estimates:
mean in group MUT  mean in group WT 
         1.954045          2.052696 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = -0.55686, df = 28.757, p-value = 0.5819
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.6595740  0.9494609
sample estimates:
mean in group MUT  mean in group WT 
        0.2369000         0.5919565 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = -3.7918, df = 43.355, p-value = 0.0004583
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -2.1732789 -0.6644011
sample estimates:
mean in group MUT  mean in group WT 
          0.71900           2.13784 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 0.93466, df = 34.257, p-value = 0.3565
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.2382042  0.6441042
sample estimates:
mean in group MUT  mean in group WT 
          3.24755           3.04460 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = 0.63465, df = 37.688, p-value = 0.5295
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4989242  0.9544242
sample estimates:
mean in group MUT  mean in group WT 
          3.35485           3.12710 


> element_stats_nonpara[[7]]
[[1]]

	Wilcoxon rank sum exact test

data:  df_fp_col[, int] by Genotype
W = 153, p-value = 0.1418
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum exact test

data:  df_fp_si[, int] by Genotype
W = 265, p-value = 0.796
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_col[, int] by Genotype
W = 231, p-value = 0.9734
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_si[, int] by Genotype
W = 134, p-value = 0.001114
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum test with continuity correction

data:  df_ts_col[, int] by Genotype
W = 224, p-value = 0.525
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum exact test

data:  df_ts_si[, int] by Genotype
W = 224, p-value = 0.5291
alternative hypothesis: true location shift is not equal to 0
```