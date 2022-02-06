###  run gJLS2 (plink plugin via PLINK1.9)

## PLUGIN scripts:
# install.packages("Rserve")
# library(Rserve)
# Rserve()

# run Rserve in the background 

R CMD Rserve --RS-port 8221

plink --bfile ./input/chrX_5_snp \
--R run_gJLS2PLINK_Xchr.R \
--pheno ./input/Pheno.txt \
--pheno-name pheno1 \
--R-port 8221 \
--covar ./input/Pheno.txt \
--covar-name SEX,covar1,covar2,covar3 \
--out ./output/testRun


### debug codes:

plink --bfile ./input/chrX_5_snp \
--R run_gJLS2PLINK_Xchr.R \
--R-debug \
--pheno ./input/Pheno.txt \
--pheno-name pheno1 \
--R-port 8221 \
--covar ./input/Pheno.txt \
--covar-name SEX,covar1,covar2,covar3 \
--out ./output/testRun


### important: SEX must always be the first covar-name specified

source("testRun.debug.R")

 require(gJLS2)
 
  f1 <- function(s) 
       {    
      r <-  gJLS2(GENO=s, SEX=ifelse(COVAR[,1]==0, 2, COVAR[,1]), Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE, genotypic=TRUE)
      rr <- as.numeric(r[3:5])
      
      c( length(rr) , rr )
       }
       
apply( GENO , 2 , f1 )
      
      
      
      