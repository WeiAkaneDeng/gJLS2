
# Using PLINK R Plugin and R scripts 

Considering the computational and memory requirement of genome-wide data, the package offers two possible approaches for running the gJLS analyses: 1) using the stand-alone package in R environment (GUI or as a command-line program through \texttt{Rscript}), or 2) using the genetic software PLINK via an R plugin function.

In general, we recommend the users to divide the genotype data by chromosome and maybe take advantage of parallel computing when using a server with multiple cores or processors.


## PLINK and an R plugin function

We illustrate gJLS analyses using simulated phenotype, covariates, and real X-chromosome genotypes from 1000 Genomes Project via a simple yet adaptable R plugin function; the function can be easily modified to suit the user's needs. 

Before we start, we have to make sure both \package{gJLS2} and \package{Rserve} packages are installed. More details on the R plugin function can be found [here](https://www.cog-genomics.org/plink/1.9/rserve). In some cases, the default port does not work and you can specify a port number when starting \package{Rserve} and in PLINK. 

The following R plugin script is quite simple since the majority of the coding is now nested in the gJLS2 function. Notice that for X-chromosome, the genetic sex must be provided and this is done through setting the first column of the covariate file to be the sex variable. To run PLINK R plugin function, save this coding chunk into an R script and give it a name, such as 'run_gJLS2PLINK_Xchr.R'.

```
Rplink <- function(PHENO,GENO,CLUSTER,COVAR){
 require(gJLS2)
 
  f1 <- function(s) 
       {    
      r <-  gJLS2(GENO=s, SEX=COVAR[,1], Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE, genotypic=TRUE)
      rr <- as.numeric(r[3:5])
      
      c( length(rr) , rr )
       }
      apply( GENO , 2 , f1 )

}
```

Notice that the first variable in the PLINK *COVAR* term is in fact the sex variable, taken from the .fam file. So please be sure to check if sex information is correctly stored in .fam file prior to this analysis.

The remaining step is done in the bash command line by calling PLINK and you can test the script with the files extracted from the file [here](https://github.com/WeiAkaneDeng/gJLS2/blob/main/inst/extdata/input.zip).


```
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

The PLINK has an option to debug the script when necessary, it can be done on a subset of the data using the following:


```bash
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
This generates a *testRun.debug.R* file in which the data needed for the analysis are stored, and the user can proceed to debug in R:


```
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

This option has the flexibility to perform data filtering via PLINK, but sometimes setting up *Rserve* to communicate with PLINK properly could be a problem on some servers. We recommend running this analysis by either chromosome or smaller units when sample sizes become large as the R *apply()* function supported by PLINK is not the most efficient. There is a risk of losing the results if the analyese were not fully completed on *GENO* due to external factors, such as power loss or server downtime.

To have better write control over larger analyses, we recommend using the Rscript option below.


## Rscript 

The Rscript [run_gJLS2.R](https://github.com/WeiAkaneDeng/gJLS2/blob/main/inst/extdata/run_gJLS2.R) in the extdata folder can serve as a starting point to customized the analyses for each user. The arguments available in this Rscript included the very basic ones and additional ones can be added easily. A useful option is "write", where the user can specify the chunk size for results to be written in the file as the program runs. Another important feature is the "-nThreads" option, where the user can take advantage of the multiple cores and processors available from high performance computing clusters.

As an example, in the extdata directory, the Rscript can be executed on the command line and output the results:


```bash
Rscript run_gJLS2.R --bfile ./input/chrX_5_snp \
--pfile ./input/Pheno.txt \
--pheno pheno1 \
--Xchr TRUE \
--nThreads 2 \
--covar SEX,covar1,covar2,covar3 \
--out ./output/testRun.results.txt
```

Alternatively, summary statistics can also be analyzed using the command-line option:


```bash
Rscript run_gJLS2.R --sumfile ./input/GIANT_BMI_chr16_gJLS_summary.txt \
--out ./output/GIANT_BMI_Sum.chr16_results.txt &
```

