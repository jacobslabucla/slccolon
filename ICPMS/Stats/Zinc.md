```R
> # Zinc
> element_stats_para[[4]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -2.4015, df = 35.745, p-value = 0.02166
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -117.211115   -9.865528
sample estimates:
mean in group MUT  mean in group WT 
         294.5423          358.0806 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = -0.44508, df = 42.996, p-value = 0.6585
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -71.91412  45.91030
sample estimates:
mean in group MUT  mean in group WT 
          207.605           220.607 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = 1.2986, df = 19.216, p-value = 0.2094
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -473.7568 2025.8584
sample estimates:
mean in group MUT  mean in group WT 
        1006.1055          230.0547 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = 0.70932, df = 41.175, p-value = 0.4821
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -88.82218 185.01275
sample estimates:
mean in group MUT  mean in group WT 
         241.6401          193.5448 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 1.2992, df = 36.659, p-value = 0.202
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -7.293193 33.338393
sample estimates:
mean in group MUT  mean in group WT 
         164.1181          151.0954 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = 0.90883, df = 37.863, p-value = 0.3692
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -12.04897  31.67697
sample estimates:
mean in group MUT  mean in group WT 
         107.4965           97.6825 


> element_stats_nonpara[[4]]
[[1]]

	Wilcoxon rank sum exact test

data:  df_fp_col[, int] by Genotype
W = 123, p-value = 0.02291
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum exact test

data:  df_fp_si[, int] by Genotype
W = 223, p-value = 0.5067
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum exact test

data:  df_muc_col[, int] by Genotype
W = 306, p-value = 0.06563
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum exact test

data:  df_muc_si[, int] by Genotype
W = 289, p-value = 0.9837
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum exact test

data:  df_ts_col[, int] by Genotype
W = 250, p-value = 0.1826
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum exact test

data:  df_ts_si[, int] by Genotype
W = 232, p-value = 0.3983
alternative hypothesis: true location shift is not equal to 0
```