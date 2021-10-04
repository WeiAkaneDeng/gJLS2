#' Levene's regression tests for variance homogeneity by SNP genotype (X-chromosome specific)
#'
#' This function takes as input the genotype of a SNP (\code{geno_one}), the genetic sex (\code{SEX}), a quantitative trait (\code{Y}) in a sample population, and possibly additional covariates, such as principal components. The function returns the scale association \emph{p}-values for each X-chromosome SNP using the generalized Levene's test designed for X-chromosome biallelic markers.
#'
#' @param geno_one the genotype of a biallelic SNP, must be a vector of 0, 1, 2's coded for the number of reference allele. Alternatively, for imputed genotypes, it could be a matrix/vector of dosage values, numerically between 0 and 2. The length/dimension of \code{geno_one} should match that of \code{Y}, and/or \code{SEX} and \code{COVAR}.
#' @param Y a vector of quantitative traits, such as human height.
#' @param COVAR optional: a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.
#' @param SEX optional: the genetic sex of individuals in the sample population, must be a vector of 1 and 2 following the default sex code is 1 for males and 2 for females in PLINK.
#' @param genotypic optional: a logical indicating whether the variance homogeneity should be tested with respect to an additively (linearly) coded or non-additively coded \code{geno_one}. The former has one less degree of freedom than the latter and is the default option. For dosage genotypes without genotypic probabilities, \code{genotypic} is forced to be \code{FALSE}.
#' @param loc_alg a character indicating the type of algorithm to compute the centre in stage 1; the value is either "OLS", corresponding to an ordinary linear regression under Gaussian assumptions to compute the mean, or "LAD", corresponding to a quantile regression to compute the median. The recommended default option is "LAD". For the quantile regression, the function calls \code{quantreg::rq} and the median is estimated using either the "fn" (smaller samples) or "sfn" (larger samples and sparse problems) algorithm depending the sample size, for more details see \code{?quantreg::rq}.
#' @param transformed a logical indicating whether the quantitative response \code{Y} should be transformed using a rank-based method to resemble a normal distribution; recommended for traits with non-symmetric distribution. The default option is \code{TRUE}.
#'
#'
#' @import stats
#' @import quantreg
#' @import methods
#'
#' @return the Levene's test regression p-value according to the model specified.
#' @export leveneRegX_per_SNP
#' @note We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
#' the variance values tend to be proportional to mean values when stratified by \code{geno_one}.
#'
#' @examples
#' N <- 1000
#' sex <- rbinom(N, 1, 0.5)+1
#' Y <- rnorm(N)
#' genDAT <- NA
#' genDAT[sex==2] <- rbinom(sum(sex==2), 2, 0.3)
#' table(genDAT, sex)
#' genDAT[sex==1] <- rbinom(sum(sex==1), 1, 0.3)
#' table(genDAT, sex)
#'
#' leveneRegX_per_SNP(geno_one=genDAT, SEX=sex, Y=Y)
#' leveneRegX_per_SNP(geno_one=genDAT, SEX=sex, Y=Y, genotypic=TRUE)
#' leveneRegX_per_SNP(geno_one=genDAT, SEX=sex, Y=Y, loc_alg="OLS")
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#' @references Deng WQ, Mao S, Kalnapenkis A, Esko T, Magi R, Pare G, Sun L. (2019) Analytical strategies to include the X-chromosome in variance heterogeneity analyses: Evidence for trait-specific polygenic variance structure. \emph{Genet Epidemiol}. \strong{43}(7):815-830. \doi{10.1002/gepi.22247}. PMID:31332826.
#' @references Gastwirth JL, Gel YR, Miao W. (2009). The Impact of Levene's Test of Equality of Variances on Statistical Theory and Practice. \emph{Statistical Science}. \strong{24}(3) 343-360, \doi{10.1214/09-STS301}.



leveneRegX_per_SNP <- function(geno_one, SEX, Y, COVAR = NULL, genotypic = FALSE, transformed=TRUE, loc_alg = "LAD"){

  ## check minimal inputs: geno_one, Y, transformed

  if (missing(geno_one))
    stop("The geno_onetype input is missing.")
  geno_one[geno_one==-9] <- NA

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (missing(SEX))
    stop("Sex information is required for X-chromosome association.")

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


  ## change SEX values
  if (sum(SEX %in% c(1,2,NA))!=length(SEX))
    stop("Please check the SEX variable, only 1, 2, and NA are plausible values.")
  SEX[SEX==2] = 0

  ## transform before adjust for covariates

  if (class(Y)!="numeric")
    stop("Please make sure the quantitaitve trait is a numeric vector.")

  if (transformed){
    Y <- inver_norm(Y)
  }

  N = length(Y)

  if (N > 500) {
    use_method = "sfn"
  } else {
    use_method = "fn"
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


  if (length(unique(geno_one_use)) < 4 & sum(geno_one_use == 1) > 1){
    imputed = FALSE
  } else {
    imputed = TRUE
  }



  if (length(table(geno_one_use))==1) {
    warning("Monomorphic SNP detected, results will be set to N/A");
    pval <- NA

  } else {

    ### stage 1 use additive models:

    VAR <- tapply(Y, SEX, sd, na.rm=TRUE)
    w <- ifelse(SEX == as.numeric(names(VAR)[1]), VAR[1], VAR[2])

    ### obtain residuals for LAD or OLS
    if (loc_alg == "LAD"){

        if (is.null(COVAR)){
        ### no covariates:
        use_dat <- complete.cases(geno_one_use, SEX, w, Y)
        datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "sex" = SEX[use_dat], "w" = w[use_dat])
        res_geno  <- try(as.numeric(resid(quantreg::rq(y~g+sex+g:sex, data = datafr, na.action=na.exclude, method=use_method))));
        datafr$d <- abs(res_geno)
        datafr$dw <- datafr$d/datafr$w

      } else {
        use_dat <- complete.cases(geno_one_use, SEX, w, Y, COVAR)
        if (is.null(dim(COVAR))){
          datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat], "sex" = SEX[use_dat], "w" = w[use_dat])
        } else {
          datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat,], "sex" = SEX[use_dat], "w" = w[use_dat])
        }
        res_geno  <- try(as.numeric(resid(quantreg::rq(as.formula(paste("y ~ g+sex+g:sex+", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data =  datafr,na.action=na.exclude, method=use_method))));
        datafr$d <- abs(res_geno)
        datafr$dw <- datafr$d/datafr$w
      }
    } else {
      if (is.null(COVAR)){
        ### no covariates:
        use_dat <- complete.cases(geno_one_use, SEX, w, Y)
        datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "sex" = SEX[use_dat], "w" = w[use_dat])
        res_geno  <- try(as.numeric(resid(lm(y~g+sex+g:sex, data = datafr, na.action=na.exclude))));
        datafr$d <- abs(res_geno)
        datafr$dw <- datafr$d/datafr$w
      } else {
        use_dat <- complete.cases(geno_one_use, SEX, w, Y, COVAR)
        if (is.null(dim(COVAR))){
          datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat], "sex" = SEX[use_dat], "w" = w[use_dat])
        } else {
          datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat,], "sex" = SEX[use_dat], "w" = w[use_dat])
        }
        res_geno  <- try(as.numeric(resid(lm(as.formula(paste("y ~ g+sex+g:sex +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data =  datafr,na.action=na.exclude))));
        datafr$d <- abs(res_geno)
        datafr$dw <- datafr$d/datafr$w
      }
    }


    #### step 2:

    if (methods::is(res_geno, "try-error")){

      pval <- NA

      } else{

        ### no dosage, genotype hard call only
        if (!imputed) {
            gL_model_null <- lm(dw~sex, data = datafr, na.action=na.exclude)
            if (genotypic){
              datafr$g1 <- ifelse(datafr$g==1, 1, 0)
              datafr$g2 <- ifelse(datafr$g==2, 1, 0)
              gL_model <- lm(dw~sex+g1+g2+g1:sex, data = datafr,na.action=na.exclude)
              pval <-	anova(gL_model, gL_model_null)$Pr[2]
            } else {
              gL_model <- lm(dw~sex+g+g:sex, data = datafr, na.action=na.exclude)
              pval <-	anova(gL_model, gL_model_null)$Pr[2]
            }
          ### dosage, genotype probabilities
        } else {

            gL_model_null <- lm(dw~sex, data = datafr,na.action=na.exclude)
            if (genotypic & !is.null(dim(geno_one))){
              use_dat <- complete.cases(datafr, geno_one)
              datafr$g1 <- geno_one[use_dat,2]
              datafr$g2 <- geno_one[use_dat,3]
              gL_model <- lm(dw~sex+g1+g2+sex:g1, data = datafr, na.action=na.exclude)
              pval <-	anova(gL_model, gL_model_null)$Pr[2]
            } else {
              gL_model <- lm(dw~sex+g+g:sex, data = datafr,na.action=na.exclude)
              pval <-	anova(gL_model, gL_model_null)$Pr[2]
            }
          }

      }

      }

  return(pval)
}
