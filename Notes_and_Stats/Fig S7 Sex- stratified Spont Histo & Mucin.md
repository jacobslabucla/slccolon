## 5 Month

Full

```R
> wilcox.test(Score~Genotype, histology)

	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 106, p-value = 0.9439
alternative hypothesis: true location shift is not equal to 0

```

Female
```R
> wilcox.test(Score~Genotype, f_df_mut)

	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 13, p-value = 0.9168
alternative hypothesis: true location shift is not equal to 0
```
Male
```R
> wilcox.test(Score~Genotype, m_df_mut)

	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 44, p-value = 0.7803
alternative hypothesis: true location shift is not equal to 0
```

## 10 month

Full 
```R
	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 30.5, p-value = 0.04992
alternative hypothesis: true location shift is not equal to 0
```

Female
```R
	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 6.5, p-value = 0.4185
alternative hypothesis: true location shift is not equal to 0
```

Male
```R
> wilcox.test(Score~Genotype, m_df_mut)

	Wilcoxon rank sum test with continuity correction

data:  Score by Genotype
W = 8.5, p-value = 0.1021
alternative hypothesis: true location shift is not equal to 0
```