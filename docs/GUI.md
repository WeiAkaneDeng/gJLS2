
# Basic options and usage in GUI

We illustrate the use of gJLS on simulated phenotype and genotypes of X-chromosomal SNPs from the 1000 Genomes Project. The data can be loaded directly:


```r
data("chrXdat")
head(chrXdat)
#>    FID     IID PAT MAT SEX  PHENOTYPE rs5983012_A rs4119090_G rs5911042_T
#> 1 1328 NA06984   0   0   1  0.6665371           0           0           2
#> 2 1328 NA06989   0   0   2 -0.2565215           0           0           0
#> 3 1330 NA12340   0   0   1  0.2901142           0           2           0
#> 4 1330 NA12341   0   0   2 -0.4528928           0           2           0
#> 5 1330 NA12342   0   0   1 -2.2885722           0           2           2
#> 6 1334 NA12144   0   0   1 -0.2326356           0           2           2
#>   rs986810_C rs180495_G
#> 1          0          0
#> 2          1          1
#> 3          2          2
#> 4          1          0
#> 5          2          0
#> 6          0          0
```

## Phenotype data

The phenotype data included in the dataset were simulated from a standard normal distribution without influence from the genetic sex, nor genotype variables. Though it is quite common to have phenotype with sex-specific distributions, which could lead to incorrect inference if ignored.

Thus, it is always a good idea to perform some exploratory data analysis prior to running *gJLS2*, to be aware of the potential non-normality trends (e.g. skewness or modality) as well as to spot any abnormalities in the distributions (e.g. any duplicates), either alone or when stratified by genetic sex, or possibly other suspected confounding factors. 


```r
library(ggplot2)
ggplot(data = chrXdat, aes(x=PHENOTYPE, fill=as.factor(SEX))) + 
  geom_histogram(aes(y = ..density..),alpha=0.2, position="identity", 
                 binwidth = 0.2) + geom_density(alpha=0.2, size=0.2) +
  ggtitle("Historagm of Quantitative Trait") + xlab("Phenotype")
```

![plot of chunk plot](figure/plot-1.png)


```r
print(summary(chrXdat$PHENOTYPE))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -2.780851 -0.657590 -0.075777 -0.006869  0.730500  2.614816
mean(chrXdat$PHENOTYPE); sd(chrXdat$PHENOTYPE)
#> [1] -0.006869224
#> [1] 0.9743726
library(moments)
skewness(chrXdat$PHENOTYPE)
#> [1] -0.004267257
kurtosis(chrXdat$PHENOTYPE)
#> [1] 2.71659
```


## Genotype data

The samples are restricted to the unrelated European subset ($n = 471$ with non-ambiguous sex information). To cover a range of minor allele frequencies (MAF), we hand-picked the following SNPs: rs5983012 (A/G), rs986810 (C/T), rs180495 (G/A), rs5911042 (T/C), and rs4119090 (G/A) that are outside of the pseudo-autosomal region (MAF calculated in females and rounded to the nearest digit, see table below).


```
#>   CHR       SNP A1    MAF
#> 1  23 rs5983012  A 0.1016
#> 2  23  rs986810  C 0.2033
#> 3  23  rs180495  G 0.3008
#> 4  23 rs5911042  T 0.4024
#> 5  23 rs4119090  G 0.4451
```

## gJLS2 in action

The location analysis automatically returns *p*-values from the recommended association model:

```r
locReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE);
#>   CHR         SNP        gL
#> 1   X rs5983012_A 0.9877674
#> 2   X rs4119090_G 0.9201569
#> 3   X rs5911042_T 0.3898029
#> 4   X  rs986810_C 0.4619165
#> 5   X  rs180495_G 0.8767590
```

The scale analysis automatically returns *p*-values from the recommended association model and assumes the residuals at the mean-stage are calculated using Least Absolute Deviation (LAD) rather than Ordinary Least Squares (OLS):


```r
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE)
#>   CHR         SNP        gS
#> 1   X rs5983012_A 0.1391909
#> 2   X rs4119090_G 0.9828430
#> 3   X rs5911042_T 0.1487017
#> 4   X  rs986810_C 0.9563390
#> 5   X  rs180495_G 0.3476929
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE, loc_alg="OLS")
#>   CHR         SNP        gS
#> 1   X rs5983012_A 0.1739062
#> 2   X rs4119090_G 0.9999999
#> 3   X rs5911042_T 0.1163023
#> 4   X  rs986810_C 0.9581589
#> 5   X  rs180495_G 0.3619056
```

The joint-location-scale analysis is then straightforward by combining the sets of **gL** and **gS** *p*-values.

```r
gJLS2(GENO=chrXdat[,7:11], Y=chrXdat$PHENOTYPE, SEX=chrXdat$SEX, Xchr=TRUE)
#>   CHR         SNP        gL        gS      gJLS
#> 1   X rs5983012_A 0.9943538 0.1198837 0.3727472
#> 2   X rs4119090_G 0.8881506 0.9794193 0.9911401
#> 3   X rs5911042_T 0.3488576 0.1514217 0.2081701
#> 4   X  rs986810_C 0.4898597 0.9244773 0.8116064
#> 5   X  rs180495_G 0.8702304 0.3619588 0.6788681
```



# Analytical scenarios and options

In this section, we explore some common analytical scenarios and illustrate the different options available in *gJLS2*, specifically on how to deal with dosage genotypes, related samples and analysis of X-chromosome SNPs. 

As a general comment, caution should be exercised to combine *p*-values using Fisher’s method when the trait under analyses were in different scales (e.g. log-transformed for location and un-transformed for scale analyses). In this case, we recommend a rank-based inverse normal transformation (*transformed*) the quantitative trait. This is also the default option in *scaleReg* for the scale association analyses and the joint analyses.

## Imputed SNPs

To improve coverage and boost power of detectection, imputation is now routinely used in genome-wide association studies. The process of imputation has been made easy by the publicly available reference panels (1000 Genomes Project) and softwares such as Impute2 (Howie et al., 2009; PLoS Genetics) and Beagle (Browning and Browning, 2009; American Journal of Human Genetics).

For a sample with $n$ individuals, the imputed genotypes can be a matrix of posterior genotype probabilities of dimension $n\times 3$ with each column corresponding to a particular genotype (denoted by $\eta(A/A)$, $\eta(A/B)$ or $\eta(B/B)$). 

Alternatively, it is possible to have dosage data that is a vector of $n$ entries whereby the dosage value for each individual is calculated as
$$
0\times \eta(A/A) + 1\times\eta(A/B) + 2\times\eta(B/B),
$$
where $\eta(A/A) + \eta(A/B) +\eta(B/B) = 1$ and 0, 1, 2 denote the additive genotype values coded for the number of alternative alleles (usually the minor allele). 

Finally, the dosage values can be converted to discrete genotypes by setting a pre-defined threshold on the confidence of the posterior probability. For example, PLINK defines the  hardcall threshold as the distance from the nearest hardcall:
$$
0.5 \times \sum_i |g_i - round(g_i)|
$$
where $g_i$ are allele dosages ranging between 0 and 2. A hardcall threshold of 0.1 (i.e. less than) is usually used to retain SNPs with better imputation quality and a hardcall threshold of 0 sets all uncertain genotypes to be missing.

In our analyses, all three types of imputed data are accepted. The program will automatically detect the status of imputation based on 1) whether the genotype data is supplied as a vector or a *matrix*/*data.frame*; 2) whether the number of distinct genotype values exceeds 4 (0, 1, 2, and/or a missing code, usually -9 or NA). 

Following Acar and Sun (2013), here we simulate the imputed data by introducing uncertainty using a Dirichlet distribution. A parameter $a \in [0,1]$ is used to for the correct genotype category and (1 − a)/2 for the other two, where $a = 1$ corresponds to a generative model with no genotype uncertainty and $a = 0.5$ corresponds to roughly 50\% of the "best-guess" genotypes will match the correct genotype groups.



```
library("MCMCpack")
#> Loading required package: coda
#> Loading required package: MASS
#> ##
#> ## Markov Chain Monte Carlo Package (MCMCpack)
#> ## Copyright (C) 2003-2021 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park
#> ##
#> ## Support provided by the U.S. National Science Foundation
#> ## (Grants SES-0350646 and SES-0350613)
#> ##
N <- 300 
geno <- rbinom(N, 2, 0.3)
a <- 0.3 ## uncertainty
genPP <- rbind(rdirichlet(sum(geno==0),c(a,(1-a)/2,(1-a)/2)), 
        rdirichlet(sum(geno==1),c((1-a)/2,a,(1-a)/2)),
        rdirichlet(sum(geno==2),c((1-a)/2,(1-a)/2,a)))
head(genPP);
#>           [,1]         [,2]        [,3]
#> [1,] 0.7757769 2.218303e-01 0.002392827
#> [2,] 0.9984017 3.013217e-07 0.001597968
#> [3,] 0.8949852 3.126178e-02 0.073752994
#> [4,] 0.3116164 3.995173e-03 0.684388445
#> [5,] 0.9932653 6.172884e-03 0.000561811
#> [6,] 0.2234763 7.686409e-01 0.007882812
summary(rowSums(genPP))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>       1       1       1       1       1       1
```

The genotypic probability *matrix*/*data.frame* can be analyzed directly, but needs to be an item in a list since the program cannot distinguish the supplied *matrix*/*data.frame* is for an individual imputed SNP or a matrix of SNPs in discrete/dosage genotype values. 

To analyzed the genotypic probabilities using a 2 degree-of-freedom (df) test, the genotypic option needs to be set to TRUE, otherwise the function automatically converts the genotypic probability *matrix*/*data.frame* to dosage data and then analyzes using a 1 df test. However, we cannot perform the original Levene's test on the genotypic probabilities as supposed to discrete genotype values.


```r
sex <- rbinom(N, 1, 0.5)+1 ## using PLINK coding
y <- rnorm(N)
covar <- matrix(rnorm(N*10), ncol=10)
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, genotypic = TRUE) ## geno probabilities
#>     SNP         gL        gS       gJLS
#> 1 SNP_1 0.09348609 0.1870613 0.08824716
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar) ## geno dosage
#>     SNP         gL        gS      gJLS
#> 1 SNP_1 0.09348609 0.2358194 0.1061426
try(gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, origLev = TRUE)) ## cannot perform Levene's test
#>     SNP         gL        gS      gJLS
#> 1 SNP_1 0.09348609 0.2358194 0.1061426
```

## Related samples

The *related=TRUE* option can be used to deal with related samples, usually indicated by the shared family ID (FID) in the PLINK tt>.fam* files. A cluster assignment can be specified, for example, using the FID as a factor. If the *clust* argument is not specified, then the entire samples are treated as a single group constrained by a single correlation structure. The default option assumes a compound symmetric correlation structure whereby all pairwise samples within the same cluster/group have the same correlation. The argument *cov.structure* can also take other standard classes of correlation structures listed in *corClasses* from *R* package *nlme*. See *?corClasses*. This option currently only applies to autosomal SNPs.


```r
gJLS2(GENO=geno, SEX=sex, Y=y, COVAR=covar, related=TRUE, clust = rep(1:3, c(N/2, N/4, N/4)))
#>   SNP        gL        gS      gJLS
#> 1 SNP 0.1601864 0.1752936 0.1284001
```


## X-chromosome SNPs

Genotypes of a X-chromosome marker is sex-dependent and thus are generated separately for males and females, but assuming the allele frequency to be the same (I used 0.3):

```r
genoX <- NA
genoX[sex==2] <- rbinom(sum(sex==2), 2, 0.3)
genoX[sex==1] <- rbinom(sum(sex==1), 1, 0.3)
table(genoX, sex)
#>      sex
#> genoX   1   2
#>     0 102  68
#>     1  45  68
#>     2   0  17
```

For X-chromosome analyses, the option *Xchr* must be set to *TRUE* as the function cannot distinguish autosomal genotype ro X-chromosome genotype data. For the pseudo-autosomal regions of X-chromosome, this option can be set to *FALSE*.


```r
locReg(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)
#>   CHR SNP        gL
#> 1   X SNP 0.8516859
```

The scale and joint analysis can be performed similarly following the default options of inverse-normal transformation (*transformed=TRUE*), using least absolute devation (LAD) to estimate residuals in the first stage (*loc_alg="LAD"*), and assuming an additive model (*genotypic=FALSE*):


```r
gJLS2(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)
#>   CHR SNP        gL        gS      gJLS
#> 1   X SNP 0.8564939 0.6859791 0.8999986
```


## Additional options

For both autosome and X-chromosome scale association, it is possible to choose between a genotypic test with 2 df (or 3 for X-chromosome with $GxS$ interaction) or 1 df (or 2 for X-chromosome with $GxS$) test. This is controlled by the option *genotypic=TRUE*.

As an additional option, the sex-stratified scale association *p*-values may also be reported by specifying *origLev=TRUE*, which gives the sex-specific (original) Levene's test *p*-value. This can then be combined with sex-stratified location association p-values for a sex-stratified gJLS analysis using the *gJLS2s* function.


```r
data("BMIsum")
gJLS2s(gL = BMIsum[,1:4], gS= BMIsum[,5])
# print only the first 20 SNPs
#>     CHR        SNP       BP      gL    gS         gJLS
#> 1    16  rs1000014 24417536 4.3e-01 0.710 0.6675263538
#> 2    16  rs1000047  8138689 1.2e-01 0.550 0.2453946354
#> 3    16  rs1000077 82782756 1.9e-01 0.710 0.4051345825
#> 4    16  rs1000078 82782587 7.6e-01 0.610 0.8199846965
#> 5    16  rs1000100 59269792 3.9e-01 0.260 0.3334723738
#> 6    16  rs1000174 14471228 4.3e-01 0.014 0.0367982615
#> 7    16  rs1000193  6747102 1.1e-01 0.440 0.1949675645
#> 8    16  rs1000454 81589567 5.2e-03 0.420 0.0155644883
#> 9    16  rs1000455 81589585 2.9e-03 0.400 0.0090008289
#> 10   16  rs1000640 69905668 6.1e-02 0.100 0.0372067457
#> 11   16  rs1000686 77262198 1.3e-01 0.098 0.0683247299
#> 12   16  rs1000711 83888758 5.7e-01 0.530 0.6637128513
#> 13   16  rs1000742 55325460 8.5e-01 0.910 0.9721577025
#> 14   16  rs1001170  5082231 5.3e-01 0.160 0.2940405856
#> 15   16  rs1001171  5082289 5.9e-01 0.180 0.3443461903
#> 16   16  rs1001362 56674858 3.7e-01 0.670 0.5936535273
#> 17   16  rs1001493 72625552 7.8e-05 0.720 0.0006058151
#> 18   16  rs1001553 82487917 9.5e-01 0.970 0.9968349305
#> 19   16  rs1001554 82487790 8.6e-02 0.960 0.2884836269
#> 20   16  rs1001608 50881317 1.3e-02 0.950 0.0666171253
```



