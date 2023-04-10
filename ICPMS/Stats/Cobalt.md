```R
> # Cobalt 
> element_stats_para[[2]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -2.5617, df = 38.972, p-value = 0.0144
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.58079347 -0.06827796
sample estimates:
mean in group MUT  mean in group WT 
         1.845750          2.170286 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = 0.53567, df = 41.922, p-value = 0.595
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1659541  0.2858790
sample estimates:
mean in group MUT  mean in group WT 
        0.6011364         0.5411739 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = 2.5712, df = 22.287, p-value = 0.01732
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.04359428 0.40581007
sample estimates:
mean in group MUT  mean in group WT 
       0.28105000        0.05634783 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = 1.3515, df = 36.713, p-value = 0.1848
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.009516902  0.047610815
sample estimates:
mean in group MUT  mean in group WT 
       0.05608696        0.03704000 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 2.349, df = 25.171, p-value = 0.02697
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.005978493 0.090821507
sample estimates:
mean in group MUT  mean in group WT 
           0.1716            0.1232 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = 1.807, df = 36.236, p-value = 0.07908
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.002106772  0.036606772
sample estimates:
mean in group MUT  mean in group WT 
          0.08165           0.06440 


> element_stats_nonpara[[2]]
[[1]]

	Wilcoxon rank sum test with continuity correction

data:  df_fp_col[, int] by Genotype
W = 129.5, p-value = 0.03691
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum test with continuity correction

data:  df_fp_si[, int] by Genotype
W = 289, p-value = 0.4202
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_col[, int] by Genotype
W = 352, p-value = 0.001648
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum test with continuity correction

data:  df_muc_si[, int] by Genotype
W = 341.5, p-value = 0.2621
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum test with continuity correction

data:  df_ts_col[, int] by Genotype
W = 273, p-value = 0.04979
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum test with continuity correction

data:  df_ts_si[, int] by Genotype
W = 251, p-value = 0.1718
alternative hypothesis: true location shift is not equal to 0
```