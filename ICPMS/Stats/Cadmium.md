
```R
> # Cadmium
> element_stats_para[[5]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -2.3259, df = 38.202, p-value = 0.02543
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.090456498 -0.006276836
sample estimates:
mean in group MUT  mean in group WT 
        0.2313000         0.2796667 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = -0.97209, df = 26.605, p-value = 0.3398
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.14193725  0.05072381
sample estimates:
mean in group MUT  mean in group WT 
        0.1980455         0.2436522 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = 0.62571, df = 23.358, p-value = 0.5376
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.02008876  0.03753224
sample estimates:
mean in group MUT  mean in group WT 
      0.013200000       0.004478261 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = -0.59291, df = 45.991, p-value = 0.5561
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.04148081  0.02260429
sample estimates:
mean in group MUT  mean in group WT 
       0.03852174        0.04796000 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 2.6793, df = 23.188, p-value = 0.01334
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.00252232 0.01957768
sample estimates:
mean in group MUT  mean in group WT 
          0.01865           0.00760 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = -0.61215, df = 22.371, p-value = 0.5466
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.03003443  0.01633443
sample estimates:
mean in group MUT  mean in group WT 
          0.03600           0.04285 


> element_stats_nonpara[[5]]
[[1]]

	Wilcoxon rank sum test with continuity correction

data:  df_fp_col[, int] by Genotype
W = 132, p-value = 0.04322
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum test with continuity correction

data:  df_fp_si[, int] by Genotype
W = 245, p-value = 0.8648
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_col[, int] by Genotype
W = 232, p-value = 0.9203
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_si[, int] by Genotype
W = 254.5, p-value = 0.4734
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum test with continuity correction

data:  df_ts_col[, int] by Genotype
W = 299, p-value = 0.00728
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum test with continuity correction

data:  df_ts_si[, int] by Genotype
W = 210, p-value = 0.797
alternative hypothesis: true location shift is not equal to 0
```