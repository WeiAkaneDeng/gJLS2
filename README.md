# gJLS2 Overview

A joint analysis of location-scale (JLS) can be a powerful tool in genome-wide association studies to uncover previously overlooked markers that influence the quantitative traits through both mean and variance (Soave et al., 2015; AJHG). The gJLS2 package offers updated software support to the existing gJLS <https://github.com/dsoave/gJLS/>), which was developed to deal specifically with sample correlation and group uncertainty (Soave and Sun, 2017; Biometrics). Although the work was originally motivated by genetic association studies, application of the proposed method, however, is not limited to this type of data. The gJLS testing framework is useful for other scientific studies where location and scale parameters are both of interest. 

The current gJLS2, as a unifying software package to gJLS, can additionally handle analyses of X-chromosomes through the recently available association testing methods for location (Chen et al., 2020; *Biostatistics*) and scale (Deng et al., 2019; *Genetic Epidemiology*). The package offers a convenient **PLINK R-plugin** script and a command-line **Rscript** for non-GUI usage, see the last section of [this overview](https://cloud.r-project.org/web/packages/gJLS2/vignettes/Introduction.html). It is our hope that this software package will encourage the analyses and sharing of scale and gJLS p-values.

### Quick start


To install either
```{r setup1}
install.packages("gJLS2")
```
or
```{r setup2}
#install.packages("devtools")
devtools::install_github("WeiAkaneDeng/gJLS2")
```


To load the library in R, simply run:
```{r setup}
library(gJLS2)
```

### Phenotype data preparation

Any unformated file with individual IDs in the first columns (to be linked to genotype file) and the phenotype/covariates in the remaining columns. Make sure to keep track of the column number for covariate **SEX** if you want to included it for X-chromosome analyses, as it needs to be treated separately from other covariates for Xchr.


### Genotype data preparation

PLINK format files are most compatible, or unformatted discrete genotype values (samples in each row and SNPs in each column).


### User (input data) scenarios

The gJSL2 analyzes each SNP at a time, thus it is straightforward to divide-and-conquer irrespective of your programming choice. Here I list some possible user scenarios and our recommendation on which approach to use. The reader is also encouraged to experiment to find the best approach according to their own computing specifications. You can find details on how to set up for each option [here](https://cloud.r-project.org/web/packages/gJLS2/vignettes/Introduction.html).

But briefly, you can either run the gJLS2 on your own data, or run just the scale analysis and combine with existing summary statistics, or obtain both summary statistics for location and scale (some might find [the recent UKB study on variance heterogeneity](https://cnsgenomics.com/software/osca/#DataResource) and the [GIANT summary statistics on varaibility for height and BMI](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#Variability_in_BMI_and_Height_Summary_Statistics) to be useful), and run the gJLS2 through the **gJLSs** function.




# A detailed guide to gJLS2

A joint analysis of location-scale (JLS) can be a powerful tool in genome-wide association studies to uncover previously overlooked markers that influence the quantitative traits through both mean and variance (Soave et al., 2015; AJHG). The gJLS2 package offers updated software support to the existing gJLS <https://github.com/dsoave/gJLS/>), which was developed to deal specifically with sample correlation and group uncertainty (Soave and Sun, 2017; Biometrics). Although the work was originally motivated by genetic association studies, application of the proposed method, however, is not limited to this type of data. The gJLS testing framework is useful for other scientific studies where location and scale parameters are both of interest. 

The current gJLS2, as a unifying software package to gJLS, can additionally handle analyses of X-chromosomes through the recently available association testing methods for location (Chen et al., 2020; *Biostatistics*) and scale (Deng et al., 2019; *Genetic Epidemiology*). The package will offers convenient **PLINK R-plugin** scripts and a command-line **Rscript** for non-GUI usage.

The package can also be conveniently used in 1) a bash script as an R plugin function; and 2) Rscript that can be used from command line. For more details, see Section [A note on genome-wide analyses using R Plugin and R scripts].

To load the library in R, simply run:
```{r setup}
library(gJLS2)
```

We illustrate the use of gJLS on simulated phenotype and genotypes of X-chromosomal SNPs from the 1000 Genomes Project. The data can be loaded directly:

```{r gJLS2data}
data("chrXdat")
head(chrXdat)
```

## Phenotype data

The phenotype data included in the dataset were simulated from a standard normal distribution without influence from the genetic sex, nor genotype variables. Though it is quite common to have phenotype with sex-specific distributions, which could lead to incorrect inference if ignored.

Thus, it is always a good idea to perform some exploratory data analysis prior to running \code{gJLS2}, to be aware of the potential non-normality trends (e.g. skewness or modality) as well as to spot any abnormalities in the distributions (e.g. any duplicates), either alone or when stratified by genetic sex, or possibly other suspected confounding factors. 

```{r plot, fig.width=8, fig.height=6}
library(ggplot2)
ggplot(data = chrXdat, aes(x=PHENOTYPE, fill=as.factor(SEX))) + 
  geom_histogram(aes(y = ..density..),alpha=0.2, position="identity", 
                 binwidth = 0.2) + geom_density(alpha=0.2, size=0.2) +
  ggtitle("Historagm of Quantitative Trait") + xlab("Phenotype")
```

```{r summary}
print(summary(chrXdat$PHENOTYPE))
mean(chrXdat$PHENOTYPE); sd(chrXdat$PHENOTYPE)
library(moments)
skewness(chrXdat$PHENOTYPE)
kurtosis(chrXdat$PHENOTYPE)
```


## Genotype data

The samples are restricted to the unrelated European subset ($n = 471$ with non-ambiguous sex information). To cover a range of minor allele frequencies (MAF), we hand-picked the following SNPs: rs5983012 (A/G), rs986810 (C/T), rs180495 (G/A), rs5911042 (T/C), and rs4119090 (G/A) that are outside of the pseudo-autosomal region (MAF calculated in females and rounded to the nearest digit, see table below).

```{r, echo=FALSE}
dt <- data.frame("CHR" = 23, "SNP" = c("rs5983012", 
                                      "rs986810",
                                      "rs180495",
                                      "rs5911042",
                                      "rs4119090"),
                 "A1" = c("A", "C", "G", "T", "G"),
                 "MAF" = c(0.1016, 0.2033, 0.3008, 0.4024, 0.4451))
print(dt)
```

## gJLS2 in action

The location analysis automatically returns *p*-values from the recommended association model:
```{r gJLS2_example_loc}
locReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE);
```

The scale analysis automatically returns *p*-values from the recommended association model and assumes the residuals at the mean-stage are calculated using Least Absolute Deviation (LAD) rather than Ordinary Least Squares (OLS):

```{r gJLS2_example_scale}
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE)
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE, loc_alg="OLS")
```

The joint-location-scale analysis is then straightforward by combining the sets of **gL** and **gS** *p*-values.
```{r gJLS2_example_gJLS}
gJLS2(GENO=chrXdat[,7:11], Y=chrXdat$PHENOTYPE, SEX=chrXdat$SEX, Xchr=TRUE)
```



# Analytical scenarios and options

In this section, we explore some common analytical scenarios and illustrate the different options available in <tt>gJLS2</tt>, specifically on how to deal with dosage genotypes, related samples and analysis of X-chromosome SNPs. 

As a general comment, caution should be exercised to combine *p*-values using Fisher’s method when the trait under analyses were in different scales (e.g. log-transformed for location and un-transformed for scale analyses). In this case, we recommend a rank-based inverse normal transformation (<tt>transformed</tt>) the quantitative trait. This is also the default option in <tt>scaleReg</tt> for the scale association analyses and the joint analyses.


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

In our analyses, all three types of imputed data are accepted. The program will automatically detect the status of imputation based on 1) whether the genotype data is supplied as a vector or a <tt>matrix</tt>/<tt>data.frame</tt>; 2) whether the number of distinct genotype values exceeds 4 (0, 1, 2, and/or a missing code, usually -9 or NA). 

Following Acar and Sun (2013), here we simulate the imputed data by introducing uncertainty using a Dirichlet distribution. A parameter $a \in [0,1]$ is used to for the correct genotype category and (1 − a)/2 for the other two, where $a = 1$ corresponds to a generative model with no genotype uncertainty and $a = 0.5$ corresponds to roughly 50\% of the "best-guess" genotypes will match the correct genotype groups.


```{r dosage}
library("MCMCpack")
N <- 300 
geno <- rbinom(N, 2, 0.3)
a <- 0.3 ## uncertainty
genPP <- rbind(rdirichlet(sum(geno==0),c(a,(1-a)/2,(1-a)/2)), 
        rdirichlet(sum(geno==1),c((1-a)/2,a,(1-a)/2)),
        rdirichlet(sum(geno==2),c((1-a)/2,(1-a)/2,a)))
head(genPP);
summary(rowSums(genPP))
```

The genotypic probability <tt>matrix</tt>/<tt>data.frame</tt> can be analyzed directly, but needs to be an item in a list since the program cannot distinguish the supplied <tt>matrix</tt>/<tt>data.frame</tt> is for an individual imputed SNP or a matrix of SNPs in discrete/dosage genotype values. 

To analyzed the genotypic probabilities using a 2 degree-of-freedom (df) test, the genotypic option needs to be set to TRUE, otherwise the function automatically converts the genotypic probability <tt>matrix</tt>/<tt>data.frame</tt> to dosage data and then analyzes using a 1 df test. However, we cannot perform the original Levene's test on the genotypic probabilities as supposed to discrete genotype values.

```{r dosage_exm1}
sex <- rbinom(N, 1, 0.5)+1 ## using PLINK coding
y <- rnorm(N)
covar <- matrix(rnorm(N*10), ncol=10)
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, genotypic = TRUE) ## geno probabilities
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar) ## geno dosage
try(gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, origLev = TRUE)) ## cannot perform Levene's test
```

## Related samples

The <tt>related=TRUE</tt> option can be used to deal with related samples, usually indicated by the shared family ID (FID) in the PLINK tt>.fam</tt> files. A cluster assignment can be specified, for example, using the FID as a factor. If the <tt>clust</tt> argument is not specified, then the entire samples are treated as a single group constrained by a single correlation structure. The default option assumes a compound symmetric correlation structure whereby all pairwise samples within the same cluster/group have the same correlation. The argument <tt>cov.structure</tt> can also take other standard classes of correlation structures listed in <tt>corClasses</tt> from <tt>R</tt> package <tt>nlme</tt>. See <tt>?corClasses</tt>. This option currently only applies to autosomal SNPs.

```{r related_samples}
gJLS2(GENO=geno, SEX=sex, Y=y, COVAR=covar, related=TRUE, clust = rep(1:3, c(N/2, N/4, N/4)))
```


## X-chromosome SNPs

Genotypes of a X-chromosome marker is sex-dependent and thus are generated separately for males and females, but assuming the allele frequency to be the same (I used 0.3):

```{r xchr}
genoX <- NA
genoX[sex==2] <- rbinom(sum(sex==2), 2, 0.3)
genoX[sex==1] <- rbinom(sum(sex==1), 1, 0.3)
table(genoX, sex)
```

For X-chromosome analyses, the option <tt>Xchr</tt> must be set to <tt>TRUE</tt> as the function cannot distinguish autosomal genotype ro X-chromosome genotype data. For the pseudo-autosomal regions of X-chromosome, this option can be set to <tt>FALSE</tt>.

```{r xchr_loc}
locReg(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)
```

The scale and joint analysis can be performed similarly following the default options of inverse-normal transformation (<tt>transformed=TRUE</tt>), using least absolute devation (LAD) to estimate residuals in the first stage (<tt>loc_alg="LAD"</tt>), and assuming an additive model (<tt>genotypic=FALSE</tt>):

```{r xchr_scale}
gJLS2(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)
```


## Additional options

For both autosome and X-chromosome scale association, it is possible to choose between a genotypic test with 2 df (or 3 for X-chromosome with $GxS$ interaction) or 1 df (or 2 for X-chromosome with $GxS$) test. This is controlled by the option <tt>genotypic=TRUE</tt>.

As an additional option, the sex-stratified scale association *p*-values may also be reported by specifying <tt>origLev=TRUE</tt>, which gives the sex-specific (original) Levene's test *p*-value. This can then be combined with sex-stratified location association p-values for a sex-stratified gJLS analysis using the <tt>gJLSs</tt> function.

```{r chr_scale}
gJLS2(GENO=geno, SEX=sex, Y=y, COVAR=covar, origLev=TRUE)
```


# A note on genome-wide analyses using R Plugin and R scripts 

Considering the computational and memory requirement of genome-wide data, the package offers two possible approaches for running the gJLS analyses: 1) using the stand-alone package in R environment (GUI or as a command-line program through \texttt{Rscript}), or 2) using the genetic software PLINK via an R plugin function.

In general, we recommend the users to divide the genotype data by chromosome and maybe take advantage of parallel computing when using a server with multiple cores or processors.

## Rscript 

The Rscript ("run_gJLS2.R") in the extdata folder can serve as a starting point to customized the analyses for each user. The arguments available in this Rscript included the very basic ones and additional ones can be added easily. A useful option is "write", where the user can specify the chunk size for results to be written in the file as the program runs. Another important feature is the "nThreads" option, where the user can take advantage of the multiple cores and processors available from high performance computing clusters.

As an example, in the extdata directory, the Rscript can be excuted on the command line and output the results:

```{bash, eval=FALSE}
Rscript run_gJLS2.R --bfile ./input/chrX_5_snp \
--pfile ./input/Pheno.txt \
--pheno pheno1 \
--Xchr TRUE \
--nThreads 2 \
--covar SEX,covar1,covar2,covar3 \
--out ./output/testRun.results.txt
```

Alternatively, summary statistics can also be analyzed using the command-line option:

```{bash, eval=FALSE}
Rscript run_gJLS2.R --sumfile ./input/GIANT_BMI_chr16_gJLS_summary.txt \
--out ./output/GIANT_BMI_Sum.chr16_results.txt &
```


## PLINK and an R plugin function

We illustrate gJLS analyses using simulated phenotype, covariates, and real X-chromosome genotypes from 1000 Genomes Project via a simple yet adaptable R plugin function; the function can be easily modified to suit the user's needs. 

Before we start, we have to make sure both \package{gJLS2} and \package{Rserve} packages are installed. More details on the R plugin function can be found [here](https://www.cog-genomics.org/plink/1.9/rserve). In some cases, the default port does not work and you can specify a port number when starting \package{Rserve} and in PLINK.

The following R plugin script is quite simple since the majority of the coding is now nested in the gJLS2 function. Notice that for X-chromosome, the genetic sex must be provided and this is done through setting the first column of the covariate file to be the sex variable. To run PLINK R plugin function, save this coding chunk into an R script and give it a name, such as 'run_gJLS2.R'.

```{r}
Rplink <- function(PHENO,GENO,CLUSTER,COVAR){
 require(gJLS2)
 
  f1 <- function(s) 
       {    
      r <-  gJLS2(GENO=s, SEX=ifelse(COVAR[,1]==0, 2, COVAR[,1]), Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE)
      rr <- as.numeric(r[3:5])
      c( length(rr) , rr )
       }
      apply( GENO , 2 , f1 )
}
```

The remaining step is done in the bash command line by calling PLINK and you can test the script with the files extracted from the file [1KG_example.zip](https://github.com/WeiAkaneDeng/Xvarhet/blob/master/1kgExample.zip).

```{bash, eval=FALSE}
R CMD Rserve --RS-port 8221 --no-save
plink --bfile ./input/chrX_5_snp \
--R run_gJLS2PLINK_Xchr.R \
--pheno ./input/Pheno.txt \
--pheno-name pheno1 \
--R-port 8221 \
--covar ./input/Pheno.txt \
--covar-name SEX,covar1,covar2,covar3 \
--out ./output/testRun
```

The PLINK only has an option to debug the script when necessary, it can be done on a subset of the data using the following:

```{bash, eval=FALSE}
plink --bfile ./input/chrX_5_snp \
--R run_gJLS2PLINK_Xchr.R \
--R-debug \
--pheno ./input/Pheno.txt \
--pheno-name pheno1 \
--R-port 8221 \
--covar ./input/Pheno.txt \
--covar-name SEX,covar1,covar2,covar3 \
--out ./output/testRun
```
This generates a "testRun.debug.R" file in which the data needed for the analysis are stored, and the user can proceed to debug in R:

```{r, eval=FALSE}
source("testRun.debug.R")
 require(gJLS2)
 
  f1 <- function(s) 
       {    
      r <-  gJLS2(GENO=s, SEX=ifelse(COVAR[,1]==0, 2, COVAR[,1]), Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE)
      rr <- as.numeric(r[3:5])
      
      c( length(rr) , rr )
       }
       
apply(GENO , 2 , f1)
```      


### User (input data) scenarios

The gJSL2 analyzes each SNP at a time, it is straightforward to divide-and-conquer irrespective of your programming choice. Here I list some possible user scenarios and our recommendation on which approach to use. The reader is also encouraged to experiment to find the best approach according to their own computing specifications. 

* *I have genotype and phenotype data and am interested in running the gJLS analyses from scratch using the R package.* 
  + *I have a relatively small dataset with < 5,000 samples and want to investigate particular markers or gene-set analysis.* <span style="color:green">The R GUI should be sufficient for this purpose. PLINK binary data can be read with "BEDMatrix" and the .raw genotypes can be read in directly as a data.frame.</span>
  + *I have a relatively small dataset with < 5,000 samples and want to perform a genome-wide analysis.* <span style="color:green">We recommend either the Rscript or PLINK R plugin as the analysis is straighforward and simple to break the job by chromosome or apply other filters within PLINK.</span>
  + *I have datasets with > 5,000 samples.* <span style="color:green">The Rscript option will work nicely on binary files using "BEDMatrix", the multiple cores option can speed up the computational time, and the user specified write size makes sure that no results are lost in the process.</span>


* *I have genotype and phenotype data, but already ran the location analysis (GWAS) using other methods (e.g. PLINK or BOLT-LMM) or plan to use published GWAS p-values. I am interested in the gS analysis and the combined gJLS analyses using existing location p-values in place of the gL.*

  + *I have genotype and phenotype data, but I am only interested in the scale analysis.*<span style="color:green">The user can run scale analysis alone using \code{scaleReg} using the preferred approach depending on sample size and computational requirements. </span>
  + *I have no individual-level data, but I have p-values for location and scale analyses from external sources (e.g. GIANT) and am interested in the gJLS analysis.* <span style="color:green">Once the gL and gS *p*-values have been obtained, the function \code{gJLS2s} can be called in R GUI to produce the combined gJLS *p*-values without access to individual-level data.</span>



