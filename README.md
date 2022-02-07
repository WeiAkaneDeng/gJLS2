
# gJLS2: A Generalized Joint Location and Scale Tool for Genetic Association Testing

A joint analysis of location-scale (JLS) can be a powerful tool in genome-wide association studies to uncover previously overlooked markers that influence the quantitative traits through both mean and variance (Soave et al., 2015; AJHG). The gJLS2 package offers updated software support to the existing [gJLS](https://github.com/dsoave/gJLS/), which was developed to deal specifically with sample correlation and group uncertainty (Soave and Sun, 2017; Biometrics). Although the work was originally motivated by genetic association studies, application of the proposed method, however, is not limited to this type of data. The gJLS testing framework is useful for other scientific studies where location and scale parameters are both of interest. 

The current gJLS2, as a unifying software package to gJLS, can additionally handle analyses of X-chromosomes through the recently available association testing methods for location (Chen et al., 2020; *Biostatistics*) and scale (Deng et al., 2019; *Genetic Epidemiology*). The package offers a convenient **PLINK R-plugin** script and a command-line **Rscript** for non-GUI usage. 

For details of the manual, please visit <https://weiakanedeng.github.io/gJLS2/>. It is our hope that this software package will encourage the analyses and sharing of scale and gJLS p-values.


## Quick start

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

## Functionalities

The package has three main functionalities 

- location analysis for X-chromosome markers (Xchr association testing)
- scale analysis for autosomal and X-chromosome markers (variance heterogeneity testing)
- joint-location-scale analysis for both individual-level data and summary statistics

Details on the methodology and recommendation for usage can be found [here](https://weiakanedeng.github.io/gJLS2/).


## Usage

The R package can be used 

- directly in GUI for smaller analysis with sample size under 5,000 and a handful of markers
- via a [PLINK R plug-in](https://www.cog-genomics.org/plink/1.9/rserve) genome-wide (per chromosome) for sample size under 5,000
- via a Rscript for larger analyses

Details on how to prepare data, set up for each options can be found [here](https://weiakanedeng.github.io/gJLS2/). 


## Contact or Reporting issues

Any issues can be reported under [Issues](https://github.com/WeiAkaneDeng/gJLS2/issues) and if not resolved within a week, please email [Wei](dengwq@mcmaster.ca).





