## DSS 
### Percent Body Weight 
Males
```R
Fixed effects:  Score ~ Genotype 
                   Value  Std.Error  DF    t-value p-value
(Intercept) -0.010684424 0.01126022 225 -0.9488646  0.3437
GenotypeHET -0.002532013 0.01584219 225 -0.1598272  0.8732
GenotypeMUT  0.025202489 0.01633011 225  1.5433146  0.1242
```
Females
```R
Fixed effects:  Score ~ Genotype 
               Value Std.Error DF   t-value p-value
(Intercept) -5.08300  1.066643 99 -4.765420  0.0000
GenotypeHET  1.59525  1.508460  8  1.057535  0.3212
GenotypeMUT  0.90100  1.629323  8  0.552990  0.5954
```

## TNBS 
### Percent Body Weight 
Males
```R
Fixed effects:  Score ~ Genotype 
                 Value Std.Error  DF    t-value p-value
(Intercept) -1.0684424  1.126022 225 -0.9488646  0.3437
GenotypeHET -0.2532013  1.584219 225 -0.1598272  0.8732
GenotypeMUT  2.5202489  1.633010 225  1.5433146  0.1242
```
Females
```R
Fixed effects:  Score ~ Genotype 
                Value Std.Error  DF    t-value p-value
(Intercept) -3.419123  1.455670 129 -2.3488317  0.0204
GenotypeHET  0.368923  2.058628  26  0.1792081  0.8592
GenotypeMUT -1.380806  2.146753  26 -0.6432069  0.5257
```
Females Day 3
```R
> pct_weight_day3 <- pct_weight_long %>% filter(Sex=="Female" & Day=="3" & Genotype!="HET")
> t.test(Score~Genotype,pct_weight_day3)

	Welch Two Sample t-test

data:  Score by Genotype
t = 0.28152, df = 12.934, p-value = 0.7828
```
Females Day 6
```R
> pct_weight_day6 <- pct_weight_long %>% filter(Sex=="Female" & Day=="6" & Genotype!="HET")
> t.test(Score~Genotype,pct_weight_day6)

	Welch Two Sample t-test

data:  Score by Genotype
t = 1.8068, df = 8.7775, p-value = 0.1051
```
Females Day 7
```R
> pct_weight_day7 <- pct_weight_long %>% filter(Sex=="Female" & Day=="7" & Genotype!="HET")
> t.test(Score~Genotype,pct_weight_day7)

	Welch Two Sample t-test

data:  Score by Genotype
t = 1.5102, df = 8.9769, p-value = 0.1654
```