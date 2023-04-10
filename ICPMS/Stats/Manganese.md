```R
> # Manganese
> element_stats_para[[6]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -3.3312, df = 37.288, p-value = 0.001959
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -50.51974 -12.31200
sample estimates:
mean in group MUT  mean in group WT 
         122.8847          154.3006 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = -0.2259, df = 38.141, p-value = 0.8225
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -46.21021  36.93156
sample estimates:
mean in group MUT  mean in group WT 
         102.6245          107.2639 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = 2.0681, df = 22.9, p-value = 0.05011
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.004922287 19.348191852
sample estimates:
mean in group MUT  mean in group WT 
        17.692200          8.020565 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = -0.28315, df = 44.755, p-value = 0.7784
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -2.480060  1.868769
sample estimates:
mean in group MUT  mean in group WT 
         6.480435          6.786080 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 2.3364, df = 29.045, p-value = 0.02658
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.4356401 6.5529599
sample estimates:
mean in group MUT  mean in group WT 
          15.9359           12.4416 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = 0.074837, df = 37.976, p-value = 0.9407
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.213984  1.307184
sample estimates:
mean in group MUT  mean in group WT 
           5.6478            5.6012 


> element_stats_nonpara[[6]]
[[1]]

	Wilcoxon rank sum exact test

data:  df_fp_col[, int] by Genotype
W = 110, p-value = 0.008417
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum exact test

data:  df_fp_si[, int] by Genotype
W = 258, p-value = 0.9194
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum exact test

data:  df_muc_col[, int] by Genotype
W = 306, p-value = 0.06563
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum exact test

data:  df_muc_si[, int] by Genotype
W = 274, p-value = 0.7904
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum exact test

data:  df_ts_col[, int] by Genotype
W = 276, p-value = 0.04018
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum exact test

data:  df_ts_si[, int] by Genotype
W = 204, p-value = 0.9254
alternative hypothesis: true location shift is not equal to 0
```