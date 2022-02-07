# gJLS2: A Generalized Joint Location and Scale Framework for Association Testing


## Overview

A joint analysis of location-scale (JLS) can be a powerful tool in genome-wide association studies to uncover previously overlooked markers that influence the quantitative traits through both mean and variance (Soave et al., 2015; *AJHG*). The gJLS2 package offers updated software support to previous versions [JLS](https://github.com/dsoave/JLS) and [gJLS](https://github.com/dsoave/gJLS), which were developed to deal specifically with sample correlation and group uncertainty (Soave and Sun, 2017; *Biometrics*). Although the work was originally motivated by genetic association studies, application of the proposed method, however, is not limited to this type of data. The JLS testing framework is useful for other scientific studies where location and scale parameters are both of interest.

The current gJLS2 can additionally handle analyses of X-chromosomes through the recently available association testing methods for location (Chen et al., 2021; *Genetic Epidemiology*) and scale (Deng et al., 2019; *Genetic Epidemiology*). The package offers a convenient PLINK R-plugin script and a command-line Rscript for non-GUI usage.


## **Citation**

For X-chromosome **location** association analysis, please cite:

- Chen et al. (2021).  *Genetic Epidemiology*. [The X factor: a robust and powerful approach to X-chromosome-inclusive whole-genome association studies.](https://pubmed.ncbi.nlm.nih.gov/34224641/)

For X-chromosome **scale** association analysis, please cite:

- Deng et al. (2019). *Genetic Epidemiology*. [Analytical strategies to include the X-chromosome in variance heterogeneity analyses: Evidence for trait-specific polygenic variance structure](https://pubmed.ncbi.nlm.nih.gov/31332826/). 

For the gJLS2 R package software, please cite:

- Deng WQ, Sun L. (2021). [gJLS2: An R package for generalized joint location and scale analysis in X-inclusive genome-wide association
studies.](https://www.biorxiv.org/content/10.1101/2021.10.11.463951v2) *bioRxiv* 2021.10.11.463951; 


## **Resources**


### Example data

The example ChrX data described in Deng WQ, Sun L. (2021) can be found in the R package by loading 

```
data("chrXdat")
head(chrXdat)
```

The PLINK binary version of the example data can be found on the github site [here](https://github.com/WeiAkaneDeng/gJLS2/blob/main/inst/extdata/input.zip). The zip file contains the PLINK binary genotype of 5 Xchr markers and simulated phenotypes; as well as summary statistics for height and bmi extracted from the GIANT consortium for chromosome 16. 


### External Summary Data 

- The GIANT consortium has a list of location (mean) and scale (variance) association summary statistics available [here](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files) for BMI and height. 


### Talks related to X-chromosome association

- IGES and Toronto Mclaughlin Centre PRS talk on [X-chromosome topic](https://github.com/LeiSunUofT/WorkshopPRS/blob/main/slides-adv-PRS-Xchr-WDeng.pdf) 




## **What's new**

+ **2022-Feb-06** User guide is live for gJLS2 at <https://weiakanedeng.github.io/gJLS2/>.

+ **2022-Jan** gJLS2 v.0.3.0 is available at <https://github.com/WeiAkaneDeng/gJLS2>.

+ **2021-Dec** gJLS2 v.0.2.0 is available from CRAN [here](https://cran.r-project.org/web/packages/gJLS2/index.html).



### **Contact** 
If you have any suggestions or questions regarding the improvement or bugs, please email Wei Q. Deng at <dengwq@mcmaster.ca>.

For general issues of usage, it is best to post under [Issues](https://github.com/WeiAkaneDeng/gJLS2/issues).

