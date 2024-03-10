### Using Maaslin 2: log-transformed, TSS-normalized models

Datasets: 
- Luminal Colon 
- Mucosal Colon 
- Cecum 
- Proximal Colon 
- Distal Colon 
### Subsets LC and MC:
##### Model 1: ASV ~ Sequencing_Run + Sex + Site + Genotype + (1/Litter) + (1/MouseID)

Luminal colon 
```R
Warning messages:
1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  unable to evaluate scaled gradient
2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
3: Model failed to converge with 1 negative eigenvalue: -1.3e+00 
```

Mucosal Colon
```R
Warning messages:
1: Model failed to converge with 1 negative eigenvalue: -2.7e-01 
2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge with max|grad| = 0.00588974 (tol = 0.002, component 1)
```

##### Model 2: ASV ~ Sequencing_Run + Sex + Litter + Site + Genotype + (1/MouseID)
Luminal colon 
```R
fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

 Mucosal colon
```R
fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

##### Model 3: ASV ~ Sequencing_Run + Sex  + Site + Genotype + (1/MouseID)

Luminal Colon
- No errors
- 0 significant features 
Mucosal Colon
- No errors
- `Enterorhabdus` increased in HET compared to WT 

##### Model 4: ASV ~ Sequencing_Run + Litter  + Site + Genotype + (1/MouseID)

Luminal Colon

```R
fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

Mucosal Colon
```R
fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

##### Model 8: ASV ~ Site + Sex + Genotype + (1/MouseID) + (1/Litter)

Luminal Colon
- HET have 19 features q <0.25, 0 features with q<0.05
	- Increased
		- Bacteroides ovatus
		- Parasutterella
		- Candidatus Arthromitus
	- Reduced 
		- Roseburia
		- Gastranaerophilales
		- Colidextribacter 
- MUT have 17 features with q < 0.25, 2 with q<0.05* Normal visualization (and pick out a couple promising guys)
	- Increased 
		- Sphingobium yanoikuyae
		- Ruminococcaeceae
		- UCG-010
		- Colidextribacter **
		- Intestinimonas
		- Roseburia
		- [[Marvinbryantia **]]
		- Oscillibacter
		- Lachnospiraceae

Mucosal Colon
- HET have 24 features q< 0.25
	- Increased
		- Enterorhabdus
		- Bacteroides caecimuris
		- Lachnospiraceae AC2044
		- Candidatus Arthromitus
		- Lachnospiraceae ASV 
	- Reduced 
		- Akkermansia muciniphila
		- Incertae Sedis
		- Butyricicoccus
		- Odoribacter
		- Muribaculaceae (2 ASV)
		- Lachnoclostridium 
		- Lachnospiraceae 
- MUT have 19 features with q < 0.25 and only 1 with q<0.05*
	- Increased 
		- Enterorhabdus
		- Lachnospiraceae ASV 
		- Incertae Sedis
		- Erysipelatoclostridium
		- Lachnospiraceae NK4A136 group
		- [[Marvinbryantia]]
		- Lachnospiraceae ASV ** 
	- Reduced
		- Romboutsia
		- Alistipes 
		- Bacillus
		- Erysipelotrichaceae 
		- Lots of Lachnospiraceae 

### Subsets Cecum, Proximal Colon, and Distal Colon
##### Model 3: ASV ~ Sequencing_Run + Sex  + SampleType + Genotype + (1/Litter)

Distal Colon
- HET have increased Bacteroides intestinalis compared to WT 
- MUT have increased Marvinbryantia compared to WT 

Cecum 
- MUT have 4 DAT q<0.05
	- increased Marvinbryantia
	- increased Lachnospiraceae ASV 
	- reduced Alistipes
	- reduced Lachnospiraceae NK4A136 group 
- HET have 5 DAT 
	- reduced Oscillospiraceae
	- reduced Ruminococcaceae
	- reduced Lachnospiraceae x2 ASVs 
	- reduced Lachnoclostridium

Proximal Colon
- MUT have 9 DAT 
	- reduced Lachnospiraceae A2 
	- reduced Lachnospiraceae 
	- increased Marvinbryantia
	- reduced Erysipelotrichaceae
	- reduced Anaerostipes
	- reduced Lachnospiraceae
	- reduced Lachnospiraceae NK4136
	- reduced Lachnospiraceae ASF356 
	- reduced Alistipes
- HET have 4 DAT 
	- reduced Incertae Sedis
	- increased Lachnospiraceae AC2044
	- reduced Lachnoclostridium
	- increased Candidatus Arthromitus

### Subsets Luminal Cecum, Luminal Proximal Colon, and Luminal Distal Colon