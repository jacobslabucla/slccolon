```R
> # Copper
> element_stats_para[[3]]
[[1]]

	Welch Two Sample t-test

data:  df_fp_col[, int] by Genotype
t = -2.6207, df = 35.498, p-value = 0.01283
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -13.920672  -1.771137
sample estimates:
mean in group MUT  mean in group WT 
          33.4770           41.3229 


[[2]]

	Welch Two Sample t-test

data:  df_fp_si[, int] by Genotype
t = 0.77418, df = 40.939, p-value = 0.4433
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -3.361766  7.541146
sample estimates:
mean in group MUT  mean in group WT 
         21.20186          19.11217 


[[3]]

	Welch Two Sample t-test

data:  df_muc_col[, int] by Genotype
t = 1.1339, df = 19.158, p-value = 0.2708
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -183.859  619.092
sample estimates:
mean in group MUT  mean in group WT 
        278.54640          60.92987 


[[4]]

	Welch Two Sample t-test

data:  df_muc_si[, int] by Genotype
t = 0.33436, df = 45.298, p-value = 0.7396
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -23.09221  32.28754
sample estimates:
mean in group MUT  mean in group WT 
         51.43530          46.83764 


[[5]]

	Welch Two Sample t-test

data:  df_ts_col[, int] by Genotype
t = 1.4586, df = 37.744, p-value = 0.153
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6520902  4.0111902
sample estimates:
mean in group MUT  mean in group WT 
         16.36765          14.68810 


[[6]]

	Welch Two Sample t-test

data:  df_ts_si[, int] by Genotype
t = 0.76343, df = 37.034, p-value = 0.45
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.9706409  2.1443409
sample estimates:
mean in group MUT  mean in group WT 
          7.83725           7.25040 


> element_stats_nonpara[[3]]
[[1]]

	Wilcoxon rank sum exact test

data:  df_fp_col[, int] by Genotype
W = 123, p-value = 0.02291
alternative hypothesis: true location shift is not equal to 0


[[2]]

	Wilcoxon rank sum exact test

data:  df_fp_si[, int] by Genotype
W = 276, p-value = 0.6128
alternative hypothesis: true location shift is not equal to 0


[[3]]

	Wilcoxon rank sum exact test

data:  df_muc_col[, int] by Genotype
W = 305, p-value = 0.06937
alternative hypothesis: true location shift is not equal to 0


[[4]]

	Wilcoxon rank sum exact test

data:  df_muc_si[, int] by Genotype
W = 314, p-value = 0.5948
alternative hypothesis: true location shift is not equal to 0


[[5]]

	Wilcoxon rank sum exact test

data:  df_ts_col[, int] by Genotype
W = 262, p-value = 0.0965
alternative hypothesis: true location shift is not equal to 0


[[6]]

	Wilcoxon rank sum exact test

data:  df_ts_si[, int] by Genotype
W = 215, p-value = 0.698
alternative hypothesis: true location shift is not equal to 0
```