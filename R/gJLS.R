#' A Generalized Joint-Location-Scale (gJLS) Test
#'
#' This function takes as input the genotype of a SNP (\code{GENO}), the SEX (\code{SEX}), and a quantitative trait (\code{Y}) in a sample population, and possibly additional covariates, such as principal components. The function returns the location and scale association \emph{p}-values for each SNP, as well as the gJLS p-value, which provides the combined evidence via Fisher's method (Soave et al., 2015, 2017). To perform this analysis genome-wide, we recommend to use the R-plugin written for PLINK, see \code{gJLSPLINK} for more details.
#'
#' @param GENO a list of a genotype matrix/vector of SNPs, must contain values 0, 1, 2's coded for the number of reference allele. Alternatively, for imputed genotypes, it could either be a vector of dosage values between 0 and 2, or a list of matrix of genotype probabilities, numerically between 0 and 1 for each genotype. The length/dimension of \code{GENO} should match that of \code{Y}, and/or \code{SEX} and \code{COVAR}.
#' @param Y a vector of quantitative traits, such as human height.
#' @param COVAR optional: a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.
#' @param SEX optional: the genetic sex of individuals in the sample population, must be a vector of 1 and 2 following the default sex code is 1 for males and 2 for females in PLINK.
#' @param genotypic a logical indicating whether the variance homogeneity should be tested with respect to an additively (linearly) coded or non-additively coded \code{geno_one}. The former has one less degree of freedom than the latter and is the default option. For dosage genotypes without genotypic probabilities, \code{genotypic} is forced to be \code{FALSE}.
#' @param loc_alg a character indicating the type of algorithm to compute the centre in stage 1; the value is either "OLS", corresponding to an ordinary linear regression under Gaussian assumptions to compute the mean, or "LAD", corresponding to a quantile regression to compute the median. The recommended default option is "LAD". For the quantile regression, the function calls \code{quantreg::rq} and the median is estimated using either the "br" (smaller samples) or "sfn" (larger samples and sparse problems) algorithm depending the sample size, for more details see \code{?quantreg::rq}.
#' @param transformed a logical indicating whether the quantitative response \code{Y} should be transformed using a rank-based method to resemble a normal distribution; recommended for traits with non-symmetric distribution. The default option is \code{TRUE}.
#' @param Xchr a logical indicator for whether the analysis is for X-chromosome SNPs.
#' @param origLev a logical indicator for whether the reported p-values should also include original Levene's test.
#' @param centre a character indicating whether the absolute deviation should be calculated with respect to "median" or "mean" in the traditional sex-specific and Fisher combined Levene's test p-values (three tests) for X-chromosome. The default value is "median". This option applies to sex-specific analysis using original Levene's test (i.e. when \code{regression}$$=$$\code{TRUE}).
#' @param related optional: a logical indicating whether the samples should be treated as related; if \code{TRUE} while no relatedness covariance information is given, it is then estimated under a \code{cov.structure} and assumes this structure among all within-group errors pertaining to the same pair/cluster if specified using \code{clust}. This option currently only applies to autosomal SNPs.
#' @param cov.structure optional: should be one of standard classes of correlation structures listed in \code{corClasses} from \pkg{R} package \pkg{nlme}. See \code{?corClasses}. The most commonly used option is \code{corCompSymm} for a compound symmetric correlation structure. This option currently only applies to autosomal SNPs.
#' @param clust optional: a factor indicating the grouping of samples; it should have at least two distinct values. It could be the family ID (FID) for family studies. This option currently only applies to autosomal SNPs.
#' @param XchrMethod an integer taking values 0 (reports all models), 1.1, 1.2, 2, 3, for the choice of X-chromosome location association testing models; for more details, see \code{\link{locReg}}.
#'
#' @importFrom methods is
#' @importFrom stats resid
#' @importFrom stats lm
#' @importFrom stats complete.cases
#' @importFrom stats na.exclude
#' @importFrom stats anova
#' @import quantreg
#' @import nlme
#'
#' @return a vector of location, scale and combined gJLS p-values for each SNP.
#' @export gJLS2
#'
#' @note For a genome-scan, we recommend to run this in PLINK via the plugin function \code{gJLSPLINK}, especially for large datasets and those with more than 20 covariates.
#' @note We highly recommend to quantile-normally transform \code{Y} for non-symmetrically distributed traits. This is typically done to avoid ‘scale-effect’ when the variance values tend to be proportional to mean values when stratified by \code{GENO}, as observed by Pare et al. (2010) and Yang et al. (2011).
#' @note For the moment, only quantitative trait \code{Y} is accepted as the subsequent generalized joint location scale (gJLS) analyses require the variance be calculated on quantitative traits. However, we are working on to include binary response for the generalized JLS analyses in the next update of gJLS.
#'
#' @examples
#' N <- 1000
#' genDAT <- rbinom(N, 2, 0.3)
#' sex <- rbinom(N, 1, 0.5)+1
#' y <- rnorm(N)
#' covar <- matrix(rnorm(N*10), ncol=10)
#'
#' gJLS2(GENO=data.frame("SNP1" = genDAT, "aSNP1" = genDAT), SEX=sex, Y=y, COVAR=covar)
#'
#' gJLS2(GENO=genDAT, SEX=sex, Y=y, COVAR=covar , Xchr=TRUE)
#'
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#' @references Soave D, Corvol H, Panjwani N, Gong J, Li W, Boëlle PY, Durie PR, Paterson AD, Rommens JM, Strug LJ, Sun L. (2015). A Joint Location-Scale Test Improves Power to Detect Associated SNPs, Gene Sets, and Pathways. \emph{American Journal of Human Genetics}. 2015 Jul 2;\strong{97}(1):125-38. \doi{10.1016/j.ajhg.2015.05.015}. PMID: 26140448; PMCID: PMC4572492.



gJLS2 <- function(GENO, Y, COVAR = NULL, SEX = NULL, Xchr = FALSE, transformed=TRUE, loc_alg = "LAD", related = FALSE, cov.structure = "corCompSymm", clust = NULL, genotypic = FALSE,  origLev = FALSE, centre = "median", XchrMethod = 3){

  if (missing(GENO))
    stop("The geno_onetype input is missing.")

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (class(Y)!="numeric")
    stop("Please make sure the quantitaitve trait is a numeric vector.")

  suppressWarnings(locP <- locReg(GENO=GENO, Y = Y, SEX=SEX, COVAR=COVAR, Xchr=Xchr, XchrMethod = XchrMethod, transformed = transformed, related = related, cov.structure = cov.structure, clust = clust))

  suppressWarnings(scalP <- scaleReg(GENO=GENO, Y = Y, COVAR = COVAR, SEX = SEX, Xchr = Xchr, transformed=transformed, loc_alg = loc_alg, related = related, cov.structure = cov.structure, clust = clust, genotypic = genotypic,  origLev = origLev, centre = centre))

  suppressMessages(merged_dat <- plyr::join(locP, scalP))

  merged_dat$gJLS <- 1-pchisq(-2*log(merged_dat$gL)-2*log(merged_dat$gS), 4)
  rownames(merged_dat) <- NULL

  return(merged_dat)
}


