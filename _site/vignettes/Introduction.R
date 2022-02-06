## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  autodep = TRUE
)

## ----setup1, eval=FALSE-------------------------------------------------------
#  install.packages("gJLS2")

## ----setup2, eval=FALSE-------------------------------------------------------
#  #install.packages("devtools")
#  devtools::install_github("WeiAkaneDeng/gJLS2")

## ----setup--------------------------------------------------------------------
library(gJLS2)

## ----gJLS2data----------------------------------------------------------------
data("chrXdat")
head(chrXdat)

## ----plot, fig.width=8, fig.height=6------------------------------------------
library(ggplot2)
ggplot(data = chrXdat, aes(x=PHENOTYPE, fill=as.factor(SEX))) + 
  geom_histogram(aes(y = ..density..),alpha=0.2, position="identity", 
                 binwidth = 0.2) + geom_density(alpha=0.2, size=0.2) +
  ggtitle("Historagm of Quantitative Trait") + xlab("Phenotype")

## ----summary------------------------------------------------------------------
print(summary(chrXdat$PHENOTYPE))
mean(chrXdat$PHENOTYPE); sd(chrXdat$PHENOTYPE)
library(moments)
skewness(chrXdat$PHENOTYPE)
kurtosis(chrXdat$PHENOTYPE)

## ---- echo=FALSE--------------------------------------------------------------
dt <- data.frame("CHR" = 23, "SNP" = c("rs5983012", 
                                      "rs986810",
                                      "rs180495",
                                      "rs5911042",
                                      "rs4119090"),
                 "A1" = c("A", "C", "G", "T", "G"),
                 "MAF" = c(0.1016, 0.2033, 0.3008, 0.4024, 0.4451))
print(dt)

## ----gJLS2_example_loc--------------------------------------------------------
locReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE);

## ----gJLS2_example_scale------------------------------------------------------
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE)
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE, loc_alg="OLS")

## ----gJLS2_example_gJLS-------------------------------------------------------
gJLS2(GENO=chrXdat[,7:11], Y=chrXdat$PHENOTYPE, SEX=chrXdat$SEX, Xchr=TRUE)

## ----dosage-------------------------------------------------------------------
library("MCMCpack")
N <- 300 
geno <- rbinom(N, 2, 0.3)
a <- 0.3 ## uncertainty
genPP <- rbind(rdirichlet(sum(geno==0),c(a,(1-a)/2,(1-a)/2)), 
        rdirichlet(sum(geno==1),c((1-a)/2,a,(1-a)/2)),
        rdirichlet(sum(geno==2),c((1-a)/2,(1-a)/2,a)))
head(genPP);
summary(rowSums(genPP))

## ----dosage_exm1--------------------------------------------------------------
sex <- rbinom(N, 1, 0.5)+1 ## using PLINK coding
y <- rnorm(N)
covar <- matrix(rnorm(N*10), ncol=10)
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, genotypic = TRUE) ## geno probabilities
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar) ## geno dosage
try(gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, origLev = TRUE)) ## cannot perform Levene's test

## ----related_samples----------------------------------------------------------
gJLS2(GENO=geno, SEX=sex, Y=y, COVAR=covar, related=TRUE, clust = rep(1:3, c(N/2, N/4, N/4)))

## ----xchr---------------------------------------------------------------------
genoX <- NA
genoX[sex==2] <- rbinom(sum(sex==2), 2, 0.3)
genoX[sex==1] <- rbinom(sum(sex==1), 1, 0.3)
table(genoX, sex)

## ----xchr_loc-----------------------------------------------------------------
locReg(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)

## ----xchr_scale---------------------------------------------------------------
gJLS2(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)

## ----chr_scale----------------------------------------------------------------
data("BMIsum")
gJLS2s(gL = BMIsum[,1:4], gS= BMIsum[,5])

## -----------------------------------------------------------------------------
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

## ---- eval=FALSE--------------------------------------------------------------
#  source("testRun.debug.R")
#  
#   require(gJLS2)
#  
#    f1 <- function(s)
#         {
#        r <-  gJLS2(GENO=s, SEX=ifelse(COVAR[,1]==0, 2, COVAR[,1]), Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE)
#        rr <- as.numeric(r[3:5])
#  
#        c( length(rr) , rr )
#         }
#  
#  apply(GENO , 2 , f1)

