## produce the correlation matrix under specification
output_correlation <- function(y, clust = NULL, cov.structure = "corCompSymm"){

  ## assuming cov.structure is either compound symmetry (cs) or ar1
  if (!cov.structure %in% c("corAR1", "corARMA", "corCAR1", "corCompSymm", "corExp", "corGaus", "corLin", "corRatio", "corSpher", "corSymm"))
    stop("The cov.structure argument does match one of the standard classes of correlation structures in corClasses. See ?nlme::corClasses for more details.")

  if (cov.structure == "corAR1") {
     correlation = nlme::corAR1(form = ~ 1 | clust)
  }

  if (cov.structure == "corARMA") {
    correlation = nlme::corARMA(form = ~ 1 | clust)
  }

  if (cov.structure == "corCAR1") {
    correlation = nlme::corCAR1(form = ~ 1 | clust)
  }

  if (cov.structure == "corCompSymm") {
    correlation = nlme::corCompSymm(form = ~ 1 | clust)
  }

  if (cov.structure == "corExp") {
    correlation = nlme::corExp(form = ~ 1 | clust)
  }

  if (cov.structure == "corGaus") {
    correlation = nlme::corGaus(form = ~ 1 | clust)
  }

  if (cov.structure == "corLin") {
    correlation = nlme::corLin(form = ~ 1 | clust)
  }

  if (cov.structure == "corRatio") {
    correlation = nlme::corRatio(form = ~ 1 | clust)
  }

  if (cov.structure == "corSpher") {
    correlation = nlme::corSpher(form = ~ 1 | clust)
  }

  if (cov.structure == "corSymm") {
     correlation = nlme::corSymm(form = ~ 1 | clust)
  }
  return(correlation)
}


## rank-based inverse normal transformation
inver_norm <- function(x){
   y <- qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))
   y <- (y-mean(y, na.rm=T)/sd(y, na.rm=T));
   y
   }
