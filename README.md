# gJLS2 Overview

A joint analysis of location-scale (JLS) can be a powerful tool in genome-wide association studies to uncover previously overlooked markers that influence the quantitative traits through both mean and variance (Soave et al., 2015; AJHG). The gJLS2 package offers updated software support to the existing gJLS <https://github.com/dsoave/gJLS/>), which was developed to deal specifically with sample correlation and group uncertainty (Soave and Sun, 2017; Biometrics). Although the work was originally motivated by genetic association studies, application of the proposed method, however, is not limited to this type of data. The gJLS testing framework is useful for other scientific studies where location and scale parameters are both of interest. 

The current gJLS2, as a unifying software package to gJLS, can additionally handle analyses of X-chromosomes through the recently available association testing methods for location (Chen et al., 2020; *Biostatistics*) and scale (Deng et al., 2019; *Genetic Epidemiology*). The package offers a convenient **PLINK R-plugin** script and a command-line **Rscript** for non-GUI usage, see the last section of [this overview](https://cloud.r-project.org/web/packages/gJLS2/vignettes/Introduction.html). It is our hope that this software package will encourage the analyses and sharing of scale and gJLS p-values.

# Quick start


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
## Overview

Please find the overview [here](https://cloud.r-project.org/web/packages/gJLS2/vignettes/Introduction.html).


## Phenotype data preparation

Any unformated file with individual IDs in the first columns (to be linked to genotype file) and the phenotype/covariates in the remaining columns. Make sure to keep track of the column number for covariate **SEX** if you want to included it for X-chromosome analyses, as it needs to be treated separately from other covariates for Xchr.


## Genotype data preparation

PLINK format files are most compatible, or unformatted discrete genotype values (samples in each row and SNPs in each column).


## User (input data) scenarios

The gJSL2 analyzes each SNP at a time, thus it is straightforward to divide-and-conquer irrespective of your programming choice. Here I list some possible user scenarios and our recommendation on which approach to use. The reader is also encouraged to experiment to find the best approach according to their own computing specifications. You can find details on how to set up for each option [here](https://cloud.r-project.org/web/packages/gJLS2/vignettes/Introduction.html).

But briefly, you can either run the gJLS2 on your own data, or run just the scale analysis and combine with existing summary statistics, or obtain both summary statistics for location and scale (some might find [the recent UKB study on variance heterogeneity](https://cnsgenomics.com/software/osca/#DataResource) and the [GIANT summary statistics on varaibility for height and BMI](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#Variability_in_BMI_and_Height_Summary_Statistics) to be useful), and run the gJLS2 through the **gJLSs** function.



