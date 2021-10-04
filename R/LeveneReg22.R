#' The generalized Levene's test via a two-stage regression for variance homogeneity by SNP genotype (autosomes)
#'
#' This function takes as input the genotype of a SNP (\code{GENO}), a quantitative trait (\code{Y}) in a sample population, and possibly additional covariates, such as principal components. The function returns the scale association \emph{p}-values for each autosomal SNP using the generalized Levene's test.
#'
#' @param geno_one the genotype of a biallelic SNP, must be a vector of 0, 1, 2's coded for the number of reference allele. Alternatively, for imputed genotypes, it could be a matrix/vector of dosage values, numerically between 0 and 2. The length/dimension of \code{geno_one} should match that of \code{Y}, and/or \code{SEX} and \code{COVAR}.
#' @param Y a vector of quantitative traits, such as human height.
#' @param COVAR optional: a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.
#' @param genotypic optional: a logical indicating whether the variance homogeneity should be tested with respect to an additively (linearly) coded or non-additively coded \code{geno_one}. The former has one less degree of freedom than the latter and is the default option. For dosage data without genotypic probabilities, \code{genotypic} is forced to be \code{FALSE}.
#' @param related optional: a logical indicating whether the samples should be treated as related; if \code{TRUE} while no relatedness covariance information is given, it is then estimated under a \code{cov.structure} and assumes this structure among all within-group errors pertaining to the same pair/cluster if specified using \code{clust}. This option currently only applies to autosomal SNPs.
#' @param cov.structure optional: should be one of standard classes of correlation structures listed in \code{corClasses} from \pkg{R} package \pkg{nlme}. See \code{?corClasses}. The most commonly used option is \code{corCompSymm} for a compound symmetric correlation structure. This option currently only applies to autosomal SNPs.
#' @param clust optional: a factor indicating the grouping of samples; it should have at least two distinct values. It could be the family ID (FID) for family studies. This option currently only applies to autosomal SNPs.
#' @param loc_alg a character indicating the type of algorithm to compute the centre in stage 1; the value is either "OLS", corresponding to an ordinary linear regression under Gaussian assumptions to compute the mean, or "LAD", corresponding to a quantile regression to compute the median. The recommended default option is "LAD". For the quantile regression, the function calls \code{quantreg::rq} and the median is estimated using either the "fn" (smaller samples) or "sfn" (larger samples and sparse problems) algorithm depending the sample size, for more details see \code{?quantreg::rq}.
#' @param transformed a logical indicating whether the quantitative response \code{Y} should be transformed using a rank-based method to resemble a normal distribution; recommended for traits with non-symmetric distribution. The default option is \code{TRUE}.
#'
#' @import stats
#' @import quantreg
#' @import methods
#' @import nlme
#'
#' @return Levene's test regression p-values for autosomal SNPs according to the model specified.
#' @export leveneRegA_per_SNP
#' @note We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
#' the variance values tend to be proportional to mean values when stratified by \code{geno_one}.
#' @note When the relatedness option is used, the computational time is expected to be longer for larger sample size ($$n > 1000$$), thus we recommend this option for smaller studies rather than large population based studies.
#' @note There is no explicit argument to supply \code{SEX} for autosomal SNPs, the user can choose to include the genetic sex of individuals as a column of the \code{COVAR} argument.
#'
#'
#' @examples
#' N <- 100
#' genDAT <- rbinom(N, 2, 0.3)
#' Y <- rnorm(N)
#' covar <- matrix(rnorm(N*10), ncol=10)
#'
#' # vanilla example:
#' leveneRegA_per_SNP(geno_one=genDAT, Y=Y, COVAR=covar)
#'
#' # relatedness samples:
#' leveneRegA_per_SNP(geno_one=genDAT, Y=Y, COVAR=covar,
#' related=TRUE)
#' leveneRegA_per_SNP(geno_one=genDAT, Y=Y, COVAR=covar,
#' related=TRUE, clust = factor(rbinom(N, 2, 0.6)))
#'
#'
#' # dosage genotypes example:
#' library("MCMCpack")
#' a <- 0.3
#' geno <- rbinom(N, 2, 0.3)
#' a <- 0.3 ## uncertainty
#' genPP <- rbind(rdirichlet(sum(geno==0),c(a,(1-a)/2,(1-a)/2)),
#'                rdirichlet(sum(geno==1),c((1-a)/2,a,(1-a)/2)),
#'                rdirichlet(sum(geno==2),c((1-a)/2,(1-a)/2,a)))
#'
#' leveneRegA_per_SNP(geno_one=genPP, Y=Y, COVAR=covar)
#' leveneRegA_per_SNP(geno_one=genPP, Y=Y, COVAR=covar,
#' genotypic=TRUE)
#'
#' # dosage and related samples:
#' leveneRegA_per_SNP(geno_one=genPP, Y=Y, COVAR=covar,
#' related=TRUE, clust = factor(rbinom(N, 1, 0.3)))
#' leveneRegA_per_SNP(geno_one=genPP, Y=Y, COVAR=covar,
#' related=TRUE, clust = factor(rbinom(N, 1, 0.3)), genotypic=TRUE)
#'
#'
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#' @references Soave D, Sun L. (2017). A generalized Levene's scale test for variance heterogeneity in the presence of sample correlation and group uncertainty. \emph{Biometrics}. \strong{73}(3):960-971. \doi{10.1111/biom.12651}. PMID: 28099998.
#'
#' @references Gastwirth JL, Gel YR, Miao W. (2009). The Impact of Levene's Test of Equality of Variances on Statistical Theory and Practice." \emph{Statistical Science}. \strong{24}(3) 343 - 360, \doi{10.1214/09-STS301}
#'




leveneRegA_per_SNP <- function(geno_one, Y, COVAR = NULL, transformed=TRUE, loc_alg = "LAD", related = FALSE, cov.structure = "corCompSymm", clust = NULL, genotypic = FALSE){

  ## check minimal inputes: geno_one, Y, transformed,
  if (missing(geno_one))
    stop("The Genotype input is missing.")
  geno_one[geno_one==-9] <- NA

  ## if genotype probability matrix
  if (!is.null(dim(geno_one))){
    geno_check <- dim(geno_one)[2]

    if (geno_check > 3 | geno_check < 2){
      stop("Please check the genotype input, the genotype probabilitie matrix should have three columns.")
    } else {
      geno_one_use <- as.numeric(as.matrix(geno_one)%*%c(0,1,2))
    }
   }

  ## if not genotype probability matrix
  if (is.null(dim(geno_one))){
    geno_check <- is.vector(geno_one)|is.numeric(geno_one)|is.factor(geno_one)
    if (!geno_check){
      stop("Please check the genotype input, it should be of class vector, numeric or factor.")
    } else {
      geno_range <- diff(range(as.numeric(geno_one), na.rm=T))
      if (geno_range>2)
        stop("Please check the genotype values, the difference between the maximum and minimum genotype value should be less than 2.")
    }
    geno_one_use <- geno_one
  }

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (class(Y)!="numeric")
    stop("Please make sure the quantitaitve trait is a numeric vector.")

  ## transform before adjust for covariates
  if (transformed){
    Y <- inver_norm(Y)
  }


  ## check covariates
  if (!is.null(COVAR)){

    if (is.null(dim(COVAR))){

      if (length(geno_one_use)!=length(COVAR)|length(geno_one_use)!=length(Y)|length(Y)!=length(COVAR))
        stop("Make sure the inputs have the same length.")

    } else {

      if (length(geno_one_use)!=dim(COVAR)[1]|length(geno_one_use)!=length(Y)|length(Y)!=dim(COVAR)[1])
        stop("Make sure the inputs have the same length.")
    }
  }

  ## transform after adjust for covariates
  if (transformed){
    Y <- inver_norm(Y)
  }

  N <- length(Y)

  if (loc_alg == "LAD") {
    if (N > 500) {
      use_method = "sfn"
    } else {
      use_method = "fn"
    }
  }

  if (length(unique(geno_one_use)) < 4 & sum(geno_one_use == 1) > 1){
    imputed = FALSE
  } else {
    imputed = TRUE
  }


### stage 1 use additive models:

    if (length(table(geno_one_use))==1) {
      warning("Monomorphic SNP detected, results will be set to N/A");
      pval <- NA
    } else {

    ### obtain residuals for LAD or OLS
    if (loc_alg == "LAD"){
       if (is.null(COVAR)){
        ### no covariates:
        use_dat <- complete.cases(geno_one_use, Y)
        datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat])
        res_geno  <- try(as.numeric(resid(quantreg::rq(y~g, data = datafr, na.action=na.exclude, method=use_method))));
        datafr$d <- abs(res_geno)
       } else {
         use_dat <- complete.cases(geno_one_use, Y, COVAR)
         if (is.null(dim(COVAR))){
           datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat])
         } else {
           datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat,])
         }
        res_geno  <- try(as.numeric(resid(quantreg::rq(as.formula(paste("y ~ g +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data =  datafr,na.action=na.exclude, method=use_method))));
        datafr$d <- abs(res_geno)
      }
       } else {
         if (is.null(COVAR)){
           ### no covariates:
           use_dat <- complete.cases(geno_one_use, Y)
           datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat])
           res_geno  <- try(as.numeric(resid(lm(y~g, data = datafr, na.action=na.exclude))));
           datafr$d <- abs(res_geno)
         } else {
           use_dat <- complete.cases(geno_one_use, Y, COVAR)
           if (is.null(dim(COVAR))){
             datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat])
           } else {
             datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat,])
           }
           res_geno  <- try(as.numeric(resid(lm(as.formula(paste("y ~ g +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data =  datafr,na.action=na.exclude))));
           datafr$d <- abs(res_geno)
         }
        }



### no dosage, genotype hard call only
if (!imputed) {
    ### analyzed now consider relatedness
    if (!related){
      ## if not related then analyze as usual
        gL_model_null <- lm(d~1, data = datafr,na.action=na.exclude)
        if (genotypic){
        datafr$gF <- as.factor(datafr$g)
        gL_model <- lm(d~gF, data = datafr,na.action=na.exclude)
        pval <-	anova(gL_model, gL_model_null)$Pr[2]
        } else {
        gL_model <- lm(d~g, data = datafr, na.action=na.exclude)
        pval <-	anova(gL_model, gL_model_null)$Pr[2]
        }
	    } else {
	      if (is.null(clust)) {
	        clust = rep(1, N)
	        warning("No cluster assignment was given, treating the samples as a single group.")
	      }
	      clust <- as.factor(clust)

	      use_dat <- complete.cases(datafr, clust)
	      datafr$clust <- clust[use_dat]
	      correlation_est = output_correlation(y = datafr$y, clust = datafr$clust, cov.structure = cov.structure)

	      if (genotypic){
	        datafr$gF <- as.factor(datafr$g)
	        fitg <- gls(d~gF, data = datafr, correlation=correlation_est, method="ML",control=lmeControl(opt = "optim"))
	        pval <-	anova(fitg,Terms=2)[1,3]
	       } else {
	        fit <- gls(y~g, data=datafr, correlation=correlation_est, method="ML",control=lmeControl(opt = "optim"))
	        pval <- anova(fit,Terms=2)[1,3]
	       }
	     }

### dosage, genotype probabilities
  } else {

    ### analyzed now consider relatedness
    if (!related){
      ## if not related then analyze as usual
      gL_model_null <- lm(d~1, data = datafr,na.action=na.exclude)

      if (genotypic & !is.null(dim(geno_one))){
        use_dat <- complete.cases(datafr, geno_one)
        datafr$g1 <- geno_one[use_dat,2]
        datafr$g2 <- geno_one[use_dat,3]
        gL_model <- lm(d~g1+g2, data = datafr, na.action=na.exclude)
        pval <-	anova(gL_model, gL_model_null)$Pr[2]
      } else {
        gL_model <- lm(d~g, data = datafr,na.action=na.exclude)
        pval <-	anova(gL_model, gL_model_null)$Pr[2]
      }

    } else {
      if (is.null(clust)) {
        clust = rep(1, N)
        warning("No cluster assignment was given, treating the samples as a single group.")
      }
      clust <- as.factor(clust)

      use_dat <- complete.cases(datafr, clust)
      datafr$clust <- clust[use_dat]
      correlation_est = output_correlation(y = datafr$y, clust = datafr$clust, cov.structure = cov.structure)

      if (genotypic){
        use_dat <- complete.cases(datafr, geno_one)
        datafr$g1 <- geno_one[use_dat,2]
        datafr$g2 <- geno_one[use_dat,3]
        fitg <- gls(d~g1+g2, data = datafr, correlation=correlation_est, method="ML",control=lmeControl(opt = "optim"))
        pval <-	anova(fitg,Terms=2:3)[1,3]
      } else {
        fit <- gls(y~g, data=datafr, correlation=correlation_est, method="ML",control=lmeControl(opt = "optim"))
        pval <- anova(fit,Terms=2)[1,3]
      }
    }



  }
}

  return(pval)
}






