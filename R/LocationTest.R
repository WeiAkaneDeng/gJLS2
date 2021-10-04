#' Location (mean-based association) test
#'
#' This function takes as input the genotype of SNPs (\code{GENO}), the SEX (\code{SEX}), and a quantitative trait (\code{Y}) in a sample population, and possibly additional covariates, such as principal components. The function returns the location association \emph{p}-values for each SNP.
#'
#' @param GENO a list of a genotype matrix/vector of SNPs, must contain values 0, 1, 2's coded for the number of reference allele. Alternatively, for imputed genotypes, it could either be a vector of dosage values between 0 and 2, or a list of matrix of genotype probabilities, numerically between 0 and 1 for each genotype. The length/dimension of \code{GENO} should match that of \code{Y}, and/or \code{SEX} and \code{COVAR}.
#' @param SEX the genetic sex of individuals in the sample population, must be a vector of 1's and 2's following PLINK default coding, where males are coded as 1 and females 2. \strong{Optional for analysis of autosomal SNPs, but required for X-chromosome.}
#' @param Y a numeric vector of quantitative trait, such as human height.
#' @param COVAR optional: a vector or a matrix of covariates, such as age or principal components.
#' @param Xchr a logical indicator for whether the analysis is for X-chromosome SNPs, if \code{TRUE} then the following association testing model is used: Y~G+G_D+S+GxS; with p-value given by comparing Y~G+S+GxS vs. Y~S (G is the additive coded genotype; G_D is an indicator for female heterozygotes).
#' @param related optional: a logical indicating whether the samples should be treated as related; if \code{TRUE} while no relatedness covariance information is given, it is then estimated under a \code{cov.structure} and assumes this structure among all within-group errors pertaining to the same pair/cluster if specified using \code{clust}. This option currently only applies to autosomal SNPs.
#' @param cov.structure optional: should be one of standard classes of correlation structures listed in \code{corClasses} from \pkg{R} package \pkg{nlme}. See \code{?corClasses}. The most commonly used option is \code{corCompSymm} for a compound symmetric correlation structure. This option currently only applies to autosomal SNPs.
#' @param clust optional: a factor indicating the grouping of samples; it should have at least two distinct values. It could be the family ID (FID) for family studies. This option currently only applies to autosomal SNPs.
#' @param transformed a logical indicating whether the quantitative response \code{Y} should be transformed using a rank-based method to resemble a normal distribution; recommended for traits with non-symmetric distribution. The default option is \code{FALSE}.
#' @param XchrMethod an integer taking values 0 (reports all models), 1.1, 1.2, 2, 3, for the choice of X-chromosome association testing models:
#' model 1,1: Y~G (females only)
#' model 1.2: Y~G (males only)
#' model 2: Y~G+S+GxSex; with p-value given by comparing Y~G+Sex+GxSex vs. Y~Sex (the additively coded G is robust to X-chromosome inactivation uncertainty). This is also the option for dosage genotypes.
#' model 3 (recommended): Y~G+G_D+S+GxSex; with p-value given by comparing Y ~ G+G_D+Sex+GxSex vs. Y ~ Sex (G_D is an indicator for female heterozygotes, this model is robust to X-chromosome inactivation uncertainty and skewed inactivation). For imputed data in the form of genotypic probabilities, the model becomes: Y ~ G1 + G2 + G1xSex + Sex, where G1 and G2 are the genotypic probabilities for the heterozygote and alternative allele homozygote.
#'
#' @import methods
#' @import stats
#' @import quantreg
#' @import nlme
#'
#' @return a vector of location association \emph{p}-values for each SNP.
#' @export locReg
#'
#'
#' @note The choice to use a rank-based inverse normal transformation is left to the user's discretion. See XXX for a discussion on the pros and cons of quantile transformation with respect to location association.
#' @note For X-chromosome markers, when the samples consist entirely of females or males, we report only results from model 1, regardless of the \code{XchrMethod} option.
#'
#' @examples
#' N <- 100
#' genDAT <- rbinom(N, 2, 0.3)
#' sex <- rbinom(N, 1, 0.5)+1
#' y <- rnorm(N)
#' COVAR <- matrix(rnorm(N*10), ncol=10)
#'
#' locReg(GENO=genDAT, SEX=sex, Y=y, COVAR=COVAR)
#'
#' # correlated example:
#' library("MASS")
#' yy <- mvrnorm(1, mu= rep(0, N), Sigma = matrix(0.3, N, N) + diag(0.7, N))
#' locReg(GENO=list("SNP1"= genDAT, "SNP2" = genDAT[sample(1:100)]),
#' SEX=sex, Y=as.numeric(yy), COVAR=COVAR, related = TRUE,
#' clust = rep(1, 100))
#'
#' # sibpair example:
#' pairedY <- mvrnorm(N/2,rep(0,2),matrix(c(1,0.2,0.2,1), 2))
#' yy <- c(pairedY[,1], pairedY[,2])
#' locReg(GENO=list("SNP1"= genDAT, "SNP2" = genDAT[sample(1:100)]),
#' SEX=sex, Y=as.numeric(yy), COVAR=COVAR, related = TRUE,
#' clust = rep(c(1:50), 2))
#'
#' # Xchr data example:
#' genDAT1 <- rep(NA, N)
#' genDAT1[sex==1] <- rbinom(sum(sex==1), 1, 0.5)
#' genDAT1[sex==2] <-rbinom(sum(sex==2), 2, 0.5)
#' locReg(GENO=genDAT1, SEX=sex, Y=y, COVAR=COVAR, Xchr=TRUE)
#'
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#'
#' @references Chen B, Craiu RV, Sun L. (2020) Bayesian model averaging for the X-chromosome inactivation dilemma in genetic association study. \emph{Biostatistics}. \strong{21}(2):319-335. \doi{10.1093/biostatistics/kxy049}. PMID: 30247537.
#' @references Chen B, Craiu RV, Strug LJ, Sun L. (2021) The X factor: A robust and powerful approach to X-chromosome-inclusive whole-genome association studies. \emph{Genetic Epidemiology}. \doi{10.1002/gepi.22422}. PMID: 34224641.
#'


locReg <- function(GENO, Y, SEX = NULL, COVAR = NULL, Xchr=FALSE, XchrMethod = 3, transformed = FALSE, related = FALSE, cov.structure = "corCompSymm", clust = NULL){

  if (missing(GENO))
    stop("The genotype input is missing.")

  if (!("matrix" %in% class(GENO) | "data.frame" %in% class(GENO) | "list" %in% class(GENO) | "vector" %in% class(GENO) | "integer" %in% class(GENO) | "numeric" %in% class(GENO))){
    stop("Please make sure the genotype data is an object of vector, matrix, data.frame for discrete genotypes or a list for dosage genotypes.")
  }

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (class(Y)!="numeric")
    stop("Please make sure the quantitaitve trait is a numeric vector.")


  if ("list" %in% class(GENO)){

    ### multiple
    output_test <- lapply(GENO, function(gg) tryCatch(locReg_per_SNP(geno_one = gg, SEX=SEX, Y=Y, COVAR = COVAR, Xchr=Xchr, XchrMethod = XchrMethod, transformed = transformed, related = related, cov.structure = cov.structure, clust = clust), error=function(e) NULL))

    output <- output_test
    change_to_NA <- unlist(lapply(output_test, is.null))

    if (sum(change_to_NA) == length(output_test)){
      output <- list(NA)[rep(1, length(output_test))]
    } else {
      NA_out <- rep(NA, length(output[[which(!change_to_NA)[1]]]))
      names(NA_out) <- names(output[[which(!change_to_NA)[1]]])
      output[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
    }

    if (is.null(names(GENO))){
      snp_name <- paste("SNP",1:length(output), sep="_")
    } else {
      snp_name <- names(GENO)
    }

    if (Xchr) {

      outputRes <- data.frame("CHR" = "X", "SNP" = snp_name, "gL" = do.call(rbind, output))

    } else {

      outputRes <- data.frame("SNP" = snp_name, "gL" = do.call(rbind, output))

    }

  }


  if ("data.frame" %in% class(GENO) | "matrix" %in% class(GENO)){

    output_test <- apply(GENO, 2, function(gg) tryCatch(locReg_per_SNP(geno_one = gg, SEX=SEX, Y=Y, COVAR = COVAR, Xchr=Xchr, XchrMethod = XchrMethod, transformed = transformed, related = related, cov.structure = cov.structure, clust = clust), error=function(e) NULL))

    output <- output_test
    change_to_NA <- sapply(output_test, is.null)

    if (sum(change_to_NA) == length(output_test)){
      output <- rep(NA, length(output_test))
    } else {
      NA_out <- rep(NA, length(output[[which(!change_to_NA)[1]]]))
      names(NA_out) <- names(output[[which(!change_to_NA)[1]]])
      output[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
    }

    if (is.null(dim(output))){
      output <- do.call(rbind, output)
    }

   if (is.null(colnames(GENO))){
      snp_name <- paste("SNP", 1:dim(output)[1], sep="_")
    } else {
      snp_name <- colnames(GENO)
    }

    if (Xchr) {

      outputRes <- data.frame("CHR" = "X", "SNP" = snp_name, "gL" = output)

    } else {

      outputRes <- data.frame("SNP" = snp_name, "gL" = output)

    }

  }

  if ("vector" %in% class(GENO) | "integer" %in% class(GENO) |  "numeric" %in% class(GENO)) {

    output <- tryCatch(locReg_per_SNP(geno_one = GENO, SEX=SEX, Y=Y, COVAR = COVAR, Xchr=Xchr, XchrMethod = XchrMethod, transformed =transformed, related = related, cov.structure = cov.structure , clust = clust),  error=function(e) NA)


    if (is.null(names(GENO))){
      snp_name <- "SNP"
    } else {
      snp_name <- names(GENO)
    }

    if (Xchr) {

      outputRes <-  data.frame("CHR" = "X", "SNP" = snp_name, "gL" = output)

    } else {

      outputRes <-  data.frame("SNP" = snp_name, "gL" = output)

    }
}

  rownames(outputRes) <- NULL
  return(outputRes)
}



##################################################################################################################################



locReg_per_SNP <- function(geno_one, Y, SEX = NULL, COVAR = NULL, Xchr=FALSE, XchrMethod = 3, transformed=FALSE, related = FALSE, cov.structure = "corCompSymm", clust = NULL){

  ## check minimal inputes: geno_one, Y, Xchr, XchrMethods, transformed
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

  N <- length(Y)
  COVAR_use <- NULL
  ## Replaced by either COVAR, or COVAR, SEX
  ## for autosome use COVAR_use, for X-chromosome use COVAR

  if (!is.null(COVAR)){

    if (is.null(dim(COVAR))){

    if (length(geno_one_use)!=length(COVAR)|length(geno_one_use)!=length(Y)|length(Y)!=length(COVAR))
      stop("Make sure the inputs have the same length.")

    } else {

      if (length(geno_one_use)!=dim(COVAR)[1]|length(geno_one_use)!=length(Y)|length(Y)!=dim(COVAR)[1])
        stop("Make sure the inputs have the same length.")
    }

    ### if not for Xchr, we can treat SEX as part of covariate
    if (!is.null(SEX) & !Xchr){

      if (sum(SEX %in% c(1,2,NA))!=length(SEX))
        stop("Please check the SEX variable, only 1, 2, and NA are plausible values.")

      SEX[SEX==2] = 0

      if (length(geno_one_use)!=length(SEX)|length(geno_one_use)!=length(Y)|length(Y)!=length(SEX))
        stop("Make sure the inputs: SEX, Y, and geno_one have the same length.")

      if (sum(SEX==1, na.rm= TRUE) == sum(!is.na(SEX))){
        warning("Only Males detected")

      } else if (sum(SEX==0, na.rm= TRUE) == sum(!is.na(SEX))) {

        warning("Only Females detected")
      }

      COVAR_use <- cbind(COVAR, SEX)
    }
    }

  ### if not for Xchr and no covariate, we can treat SEX as covariate
  if (is.null(COVAR) & !is.null(SEX) & !Xchr){

    if (sum(SEX %in% c(1,2,NA))!=length(SEX))
      stop("Please check the SEX variable, only 1, 2, and NA are plausible values.")

    SEX[SEX==2] = 0

    if (length(geno_one_use)!=length(SEX)|length(geno_one_use)!=length(Y)|length(Y)!=length(SEX))
      stop("Make sure the inputs: SEX, Y, and geno_one have the same length.")

      if (sum(SEX==1, na.rm= TRUE) == sum(!is.na(SEX))){
      warning("Only Males detected")

      } else if (sum(SEX==0, na.rm= TRUE) == sum(!is.na(SEX))) {

      warning("Only Females detected")
    } else {

      COVAR_use <- SEX
    }
    }

### begin analysis of autosomes:

if (!Xchr) {

     if (length(table(geno_one_use))==1) {
        warning("Monomorphic SNP detected, results will be set to N/A");
        pval <- NA
        names(pval) <- "p-value"

        } else {

        ## autosomal related / not related models:
          if (!related){
            ## if not related then analyze as usual

            if (is.null(COVAR_use)){
              ### no covariates:
              use_dat <- complete.cases(geno_one_use, Y)
              datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat])
              pval <- summary(lm(y ~ g, data=datafr))$coef[2,4]

             } else {
               ### covariates:
               use_dat <- complete.cases(geno_one_use, Y, COVAR_use)
               if (is.null(dim(COVAR_use))){
                 datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR_use[use_dat])
               } else {
                 datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR_use[use_dat,])
               }
               pval <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr))$coef[2,4]
             }

          } else {

            if (is.null(clust)) {
              clust = rep(1, N)
              warning("No cluster assignment was given, treating the samples as a single group.")
            }
            clust <- as.factor(clust)

            if (missing(cov.structure))
              stop("The cov.structure input is missing, should be one of the standard classes of correlation structures in corClasses. See ?nlme::corClasses for more details. The default option is corCompSymm.")

            if (is.null(COVAR_use)){
                use_dat <- complete.cases(geno_one_use, clust, Y)
                datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "clust" = clust[use_dat])
                correlation_est = output_correlation(y = datafr$y, clust = datafr$clust, cov.structure = cov.structure)
                fit <- nlme::gls(y~g, data=datafr, correlation=correlation_est, method="ML",control=lmeControl(opt = "optim"))
                pval <- anova(fit,Terms=2)[1,3]
             } else {

               use_dat <- complete.cases(geno_one_use, Y, COVAR_use, clust)
               if (is.null(dim(COVAR_use))){
                 datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR_use[use_dat], "clust" = clust[use_dat])
               } else {
                 datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR_use[use_dat,], "clust" = clust[use_dat])
               }
               correlation_est = output_correlation(y = datafr$y, clust = datafr$clust, cov.structure = cov.structure)
               fit <- nlme::gls(as.formula(paste("y ~ g +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr, correlation=correlation_est, method="ML",control=lmeControl(opt = "optim"))
               pval <- anova(fit,Terms=2)[1,3]
             }
          }

## end of autosomal analysis
        }

## Xchr not related models:
} else {

if (length(table(geno_one_use))==1) {

    warning("Monomorphic SNP detected, results will be set to N/A");
    pval <- NA
    names(pval) <- "p-value"

  } else {

  ### need to check if there is heterozygotes in discrete genotype for Xchr
  if (length(unique(geno_one_use)) <= 4 & sum(geno_one_use == 1, na.rm=T) > 1){
    imputed = FALSE
  } else {
    imputed = TRUE
  }

  if (is.null(SEX))
    stop("The sex input is missing for X-chromosome analysis.")

  if (sum(SEX %in% c(1,2,NA))!=length(SEX))
    stop("Please check the SEX variable, only 1, 2, and NA are plausible values.")

  SEX[SEX==2] = 0

    if (is.null(COVAR)){
    ## no covariates:

        use_dat <- complete.cases(geno_one_use, SEX, Y)
        datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "sex" =  SEX[use_dat])

        ### no covariates:
        if (sum(SEX==1, na.rm=TRUE) == sum(!is.na(SEX))){
          warning("Only Males detected")
          p2 <- summary(lm(y~g, data=datafr[datafr$sex==1,]))$coef[2,4]
          pval <- p2
          names(pval) <- "model 1 (males)"

        } else if (sum(SEX==0, na.rm=TRUE) == sum(!is.na(SEX))) {
          warning("Only Females detected")
          p1 <- summary(lm(y~g, data=datafr[datafr$sex==0,]))$coef[2,4]
          pval <-p1
          names(pval) <- "model 1 (females)"

        } else {

          if (!imputed){
            ## heterzygote in females only
              datafr$gD <- ifelse(datafr$g==1 & datafr$sex == 0, 1, 0)
              p1 <- summary(lm(y~g, data=datafr[datafr$sex==0,]))$coef[2,4]
              p2 <- summary(lm(y~g, data=datafr[datafr$sex==1,]))$coef[2,4]
              m1 <- lm(y~sex+g+gD+g:sex, data=datafr)
              m11 <- lm(y~sex+g+g:sex, data=datafr)
              m2 <- lm(y~sex, data=datafr)
              p_val_01 <- anova(m11, m2)[2,6]
              p_val_02 <- anova(m1, m2)[2,6]

              pval_Xchr <-c(p1, p2, p_val_01, p_val_02)
              names(pval_Xchr) <- c("model 1 (females)", "model 1 (males)", "model 2", "gL")

              if (XchrMethod == 0) {
                pval <- pval_Xchr
              } else if (XchrMethod == 1) {
                pval <- pval_Xchr[1:2]
              } else if (XchrMethod == 2) {
                pval <- pval_Xchr[3]
              } else if (XchrMethod == 3) {
                pval <- pval_Xchr[4]
              }

              } else {

              ## imputed: but in dosage values (models 1,2 only)
              if (is.null(dim(geno_one))){

                p1 <- summary(lm(y~g, data=datafr[datafr$sex==0,]))$coef[2,4]
                p2 <- summary(lm(y~g, data=datafr[datafr$sex==1,]))$coef[2,4]
                m11 <- lm(y~sex+g+g:sex, data=datafr)
                m2 <- lm(y~sex, data=datafr)
                p_val_01 <- anova(m11, m2)[2,6]

                pval_Xchr <-c(p1, p2, p_val_01)
                names(pval_Xchr) <- c("model 1 (females)", "model 1 (males)", "gL")

                if (XchrMethod == 0) {
                  pval <- pval_Xchr
                } else if (XchrMethod == 1) {
                  pval <- pval_Xchr[1:2]
                } else if (XchrMethod == 2) {
                  pval <- pval_Xchr[3]
                }

                if (XchrMethod == 3 & imputed){
                  warning("For X-chromosome analysis, dosage genotypes will be analyzed using either method 1 (females only) or method 2 (no non-additive component)")
                  warning("Returning results for method 1 and 2...")
                  pval <- pval_Xchr
                }

              } else {

                use_dat <- complete.cases(geno_one_use, SEX, Y, geno_one)
                datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "g1" = geno_one[use_dat, 2], "g2" = geno_one[use_dat, 3], "sex" =  SEX[use_dat])

                p1 <- summary(lm(y~g, data=datafr[datafr$sex==0,]))$coef[2,4]
                p2 <- summary(lm(y~g, data=datafr[datafr$sex==1,]))$coef[2,4]
                m11 <- lm(y~sex+g+g:sex, data=datafr)
                m2 <- lm(y~sex, data=datafr)
                m1 <- lm(y~sex+g1+g2+g1:sex, data=datafr)
                p_val_02 <- anova(m1, m2)[2,6]
                p_val_01 <- anova(m11, m2)[2,6]

                pval_Xchr <-c(p1, p2, p_val_01)
                names(pval_Xchr) <- c("model 1 (females)", "model 1 (males)", "model 2", "gL")

                if (XchrMethod == 0) {
                  pval <- pval_Xchr
                } else if (XchrMethod == 1) {
                  pval <- pval_Xchr[1:2]
                } else if (XchrMethod == 2) {
                  pval <- pval_Xchr[3]
                } else if (XchrMethod == 3) {
                  pval <- pval_Xchr[4]
                }

              }

            }
        }

    } else {
    ## yes covariates:
      use_dat <- complete.cases(geno_one_use, SEX, Y, COVAR)
        if (is.null(dim(COVAR))){
          datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat], "sex" = SEX[use_dat])
          } else {
        datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "covar" = COVAR[use_dat,], "sex" = SEX[use_dat])
          }

        if (sum(SEX==1,na.rm=TRUE) == sum(!is.na(SEX))){
          warning("Only Males detected")
          p2 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr))$coef[2,4]
          pval <- p2
          names(pval) <- "model 1 (males)"

        } else if (sum(SEX==0,na.rm=TRUE) == sum(!is.na(SEX))) {
          warning("Only Females detected")
          p1 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr))$coef[2,4]
          pval <-p1
          names(pval) <- "model 1 (females)"

        } else {

          if (!imputed){

           ## heterzygote in females only
          datafr$gD <- ifelse(datafr$g==1 & datafr$sex == 0, 1, 0)

            p1 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr[datafr$sex==0,])[grepl("covar",names(datafr[datafr$sex==0,]))], collapse = "+"))), data=datafr[datafr$sex==0,]))$coef[2,4]
            p2 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr[datafr$sex==1,])[grepl("covar",names(datafr[datafr$sex==1,]))], collapse = "+"))), data=datafr[datafr$sex==1,]))$coef[2,4]
            m1 <- lm(as.formula(paste("y ~ sex+g+gD+g:sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
            m11 <- lm(as.formula(paste("y ~ sex+g+g:sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
            m2 <- lm(as.formula(paste("y ~ sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
            p_val_01 <- anova(m11, m2)[2,6]
            p_val_02 <- anova(m1, m2)[2,6]

            pval_Xchr <-c(p1, p2, p_val_01, p_val_02)
            names(pval_Xchr) <- c("model 1 (females)", "model 1 (males)", "model 2", "gL")

            if (XchrMethod == 0) {
              pval <- pval_Xchr
            } else if (XchrMethod == 1) {
              pval <- pval_Xchr[1:2]
            } else if (XchrMethod == 2) {
              pval <- pval_Xchr[3]
            } else if (XchrMethod == 3) {
              pval <- pval_Xchr[4]
            }

          } else {

            ## imputed: but in dosage values (models 1,2 only)
            if (is.null(dim(geno_one))){

              p1 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr[datafr$sex==0,])[grepl("covar",names(datafr[datafr$sex==0,]))], collapse = "+"))), data=datafr[datafr$sex==0,]))$coef[2,4]
              p2 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr[datafr$sex==1,])[grepl("covar",names(datafr[datafr$sex==1,]))], collapse = "+"))), data=datafr[datafr$sex==1,]))$coef[2,4]
              m11 <- lm(as.formula(paste("y ~ sex+g+g:sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
              m2 <- lm(as.formula(paste("y ~ sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
              p_val_01 <- anova(m11, m2)[2,6]

              pval_Xchr <-c(p1, p2, p_val_01)
              names(pval_Xchr) <- c("model 1 (females)", "model 1 (males)", "gL")

              if (XchrMethod == 0) {
                pval <- pval_Xchr
              } else if (XchrMethod == 1) {
                pval <- pval_Xchr[1:2]
              } else if (XchrMethod == 2) {
                pval <- pval_Xchr[3]
              }

              if (XchrMethod == 3 & imputed){
                warning("For X-chromosome analysis, dosage genotypes will be analyzed using either method 1 (females only) or method 2 (no non-additive component)")
                warning("Returning results for method 1 and 2...")
                pval <- pval_Xchr
              }

            } else {

              use_dat <- complete.cases(geno_one_use, SEX, Y, geno_one, COVAR)
              datafr <- data.frame("y" = Y[use_dat], "g" = geno_one_use[use_dat], "g1" = geno_one[use_dat, 2], "g2" = geno_one[use_dat, 3], "sex" =  SEX[use_dat], "covar" = COVAR[use_dat])

              p1 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr[datafr$sex==0,])[grepl("covar",names(datafr[datafr$sex==0,]))], collapse = "+"))), data=datafr[datafr$sex==0,]))$coef[2,4]
              p2 <- summary(lm(as.formula(paste("y ~ g +", paste(names(datafr[datafr$sex==1,])[grepl("covar",names(datafr[datafr$sex==1,]))], collapse = "+"))), data=datafr[datafr$sex==1,]))$coef[2,4]
              m11 <- lm(as.formula(paste("y ~ sex+g+g:sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
              m2 <- lm(as.formula(paste("y ~ sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
              m1 <- lm(as.formula(paste("y ~ sex+g1+g2+g1:sex + ", paste(names(datafr)[grepl("covar",names(datafr))], collapse = "+"))), data=datafr)
              p_val_02 <- anova(m1, m2)[2,6]
              p_val_01 <- anova(m11, m2)[2,6]

              pval_Xchr <-c(p1, p2, p_val_01, p_val_02)
              names(pval_Xchr) <- c("model 1 (females)", "model 1 (males)", "model 2", "gL")

              if (XchrMethod == 0) {
                pval <- pval_Xchr
              } else if (XchrMethod == 1) {
                pval <- pval_Xchr[1:2]
              } else if (XchrMethod == 2) {
                pval <- pval_Xchr[3]
              } else if (XchrMethod == 3) {
                pval <- pval_Xchr[4]
              }

            }



          }
        }
}

}

}
## end of analysis
  return(pval)
}

