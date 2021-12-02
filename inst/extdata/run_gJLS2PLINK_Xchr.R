###  run gJLS2 (plink plugin)

Rplink <- function(PHENO,GENO,CLUSTER,COVAR){
 require(gJLS2)
 
  f1 <- function(s) 
       {    
      r <-  gJLS2(GENO=s, SEX=COVAR[,1], Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE)
      rr <- as.numeric(r[3:5])
      
      c( length(rr) , rr )
       }
       
      apply( GENO , 2 , f1 )

}



