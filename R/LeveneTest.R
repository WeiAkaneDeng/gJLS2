#' Levene's test for variance homogeneity by SNP genotypes (sex-specific p-values)
#'
#' This function takes as input the genotype of a SNP (\code{geno_one}), the genetic sex (\code{SEX}), a quantitative trait (\code{Y}) in a sample population. The function then returns the variance heterogeneity \emph{p}-values for each sex and the overall variance heterogeneity signal using Fisher's method by combining the sex-specific results.
#'
#'
#' @param geno_one the genotype of a bi-allelic SNP, must be a vector of 0, 1, 2's coded for the number of reference allele. Alternatively, for imputed genotypes, it could be a matrix/vector of dosage values, numerically between 0 and 2. The length/dimension of \code{geno_one} should match that of \code{Y}, and/or \code{SEX} and \code{COVAR}.
#' @param Y a vector of quantitative traits, such as human height.
#' @param SEX optional: the genetic sex of individuals in the sample population, must be a vector of 1 and 2 following the default sex code is 1 for males and 2 for females in PLINK.
#' @param centre a character indicating whether the absolute deviation should be calculated with respect to "median" or "mean", the default option is "median".
#' @param transformed a logical indicating whether the quantitative response \code{Y} should be transformed using a rank-based method to resemble a normal distribution; recommended for traits with non-symmetric distribution. The default option is \code{TRUE}.

#' @return a vector of Levene's test p-values according to levels specified by \code{geno_one} in each sex and the Fisher's method to combine the sex-specific Levene's test \emph{p}-values.
#'
#' @importFrom stats pchisq
#' @importFrom stats lm
#' @importFrom stats resid
#' @importFrom stats complete.cases
#' @importFrom stats anova
#'
#' @examples
#' N <- 5000
#' sex <- rbinom(N, 1, 0.5)+1
#' genDAT <- rbinom(N, 2, 0.3)
#' y <- rnorm(N);
#'
#' genDAT[sex==2] <- rbinom(sum(sex==2), 1, 0.3)
#' table(genDAT, sex)
#' leveneTests_per_SNP(geno_one=genDAT, SEX=sex, Y=y^2, transform=TRUE)
#'
#' genDAT[sex==2] <- rbinom(sum(sex==2), 1, 0.01)
#' table(genDAT, sex)
#' leveneTests_per_SNP(geno_one=genDAT, SEX=sex, Y=y^2, transform=FALSE)
#'
#' leveneTests_per_SNP(geno_one=rep(0, N), SEX=sex, Y=y^2, transform=TRUE)
#' leveneTests_per_SNP(geno_one=rep(0, N), Y=y^2, transform=TRUE)
#'
#'
#' @export leveneTests_per_SNP
#'
#' @note We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
#' the variance values tend to be proportional to mean values when stratified by \code{G}.
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#' @references  Levene H. (1960) Robust tests for equality of variances. In \emph{Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling} eds:I. Olkin, S.G. Ghurye, W. Hoeffding, W.G. Madow & H.B.Mann, pp.278-292. Stanford: Stanford University Press.
#'



leveneTests_per_SNP <- function(geno_one, SEX=NULL, Y, centre = "median", transformed=TRUE) {

  if (missing(geno_one))
    stop("The geno_onetype input is missing.")
  geno_one[geno_one==-9] <- NA

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (is.null(dim(geno_one))){
    if (length(geno_one)!=length(Y))
      stop("Make sure the inputs have the same length.")
      } else {
      stop("A vector of discrete genotypes are expected")
  }


  ## transform after adjust for covariates
  if (transformed){
    Y <- inver_norm(Y)
  }

  if (is.null(SEX)){

    if (length(table(geno_one))==1) {

      warning("Monomorphic SNP detected");
      PVAL <- data.frame(NA, NA)
      names(PVAL) <- c("Lev", "Flagged")
      rownames(PVAL) <- NULL

    } else {

      N <- length(geno_one)
      geno_one <- factor(geno_one)
      flagged <- 0

      if (min(table(geno_one)) < 30){
        warning("The minimal genotype count is less than 30, the p-value might be inflated. We recommend removing SNPs with minimal count below 30.")
        flagged <- 1
      }

      meds <- tapply(Y, geno_one, centre, na.rm = TRUE)
      resp <- abs(Y - meds[geno_one])
      Lp <- anova(lm(resp ~ factor(geno_one)))[, c(1, 4, 5)][1, 3]

      PVAL <- data.frame(Lp, as.integer(flagged))
      names(PVAL) <- c("Lev", "Flagged")
      rownames(PVAL) <- NULL
    }

  } else {

  if (length(table(geno_one))==1) {

    warning("Monomorphic SNP detected");
    PVAL <- data.frame(NA, NA, NA, NA)
    names(PVAL) <- c("LevFemale", "LevMale", "Fisher", "Flagged")
    rownames(PVAL) <- NULL

  } else {

    if (sum(SEX %in% c(1,2,NA))!=length(SEX))
      stop("Please check the SEX variable, only 1, 2, and NA are plausible values.")

    SEX[SEX==2] = 0

    N <- length(geno_one)
    geno_one <- factor(geno_one)
    group <- factor(interaction(geno_one, SEX))
    flagged <- 0

    if (min(table(group)) < 30){
      warning("The minimal genotype count in either females/males is less than 30, the p-value might be inflated. We recommend removing SNPs with minimal count below 30.")
      flagged <- 1
    }

    meds <- tapply(Y, group, centre, na.rm = TRUE)
    resp <- abs(Y - meds[group])

      if (length(table(SEX))==1) {

        ## only one sex available

        if (sum(SEX==1, na.rm=TRUE) == sum(!is.na(SEX))){
          warning("Only Males detected")
            Lp_3G_F <- NA

          if (sum(table(geno_one[SEX==1])>2)==1) {
            warning("Monomorphic SNP detected in Males");
            Lp_2G_M <- NA

          } else {
            Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(geno_one[SEX == 1])))[, c(1, 4, 5)][1, 3]
        }
        }


        if (sum(SEX==0, na.rm=TRUE) == sum(!is.na(SEX))){
          warning("Only Females detected")
            Lp_2G_M <- NA

          if (sum(table(geno_one[SEX==0])>2)==1) {
            warning("Monomorphic SNP detected in Famles");
            Lp_3G_F <- NA

            } else {
            Lp_3G_F <- anova(lm(resp[SEX == 2] ~ factor(geno_one[SEX == 2])))[, c(1, 4, 5)][1, 3]
            }
          }

      } else {

        if (sum(table(geno_one[SEX==0])>2)==1 | sum(table(geno_one[SEX==1])>2)==1) {

          if (sum(table(geno_one[SEX==0])>2)==1) {
            warning("Monomorphic SNP detected in Famles");
           Lp_3G_F <- NA
           Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(geno_one[SEX == 1])))[, c(1, 4, 5)][1, 3]

          } else if (sum(table(geno_one[SEX==1])>2)==1) {
            warning("Monomorphic SNP detected in Males");
             Lp_2G_M <- NA
             Lp_3G_F <- anova(lm(resp[SEX == 0] ~ factor(geno_one[SEX == 0])))[, c(1, 4, 5)][1, 3]
          }
        } else {
          Lp_3G_F <- anova(lm(resp[SEX == 0] ~ factor(geno_one[SEX == 0])))[, c(1, 4, 5)][1, 3]
          Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(geno_one[SEX == 1])))[, c(1, 4, 5)][1, 3]

        }
        if (is.na(Lp_3G_F)|is.na(Lp_2G_M)) {

          Lp_Fisher5 <- NA

        } else {

          Lp_Fisher5 <- pchisq(-2 * (log(Lp_3G_F) + log(Lp_2G_M)),
                             df = 4, lower.tail = FALSE)
        }
      }

    PVAL <- data.frame(Lp_3G_F, Lp_2G_M, Lp_Fisher5, as.integer(flagged))
    names(PVAL) <- c("LevFemale", "LevMale", "Fisher", "Flagged")
    rownames(PVAL) <- NULL
  }

    }
     return(PVAL)
  }
