#' Scale (variance-based association) test
#'
#' This function takes as input the genotype of SNPs (\code{GENO}), the SEX (\code{SEX}), and a quantitative trait (\code{Y}) in a sample population, and possibly additional covariates, such as principal components. The function returns the scale association \emph{p}-values for each SNP.
#'
#' @param GENO a list of a genotype matrix/vector of SNPs, must contain values 0, 1, 2's coded for the number of reference allele. Alternatively, for imputed genotypes, it could either be a vector of dosage values between 0 and 2, or a list of matrix of genotype probabilities, numerically between 0 and 1 for each genotype. The length/dimension of \code{GENO} should match that of \code{Y}, and/or \code{SEX} and \code{COVAR}.
#' @param Y a vector of quantitative traits, such as human height.
#' @param COVAR optional: a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.
#' @param SEX optional: the genetic sex of individuals in the sample population, must be a vector of 1 and 2 following the default sex code is 1 for males and 2 for females in PLINK.
#' @param genotypic a logical indicating whether the variance homogeneity should be tested with respect to an additively (linearly) coded or non-additively coded \code{geno_one}. The former has one less degree of freedom than the latter and is the default option. For dosage genotypes without genotypic probabilities, \code{genotypic} is forced to be \code{FALSE}.
#' @param loc_alg a character indicating the type of algorithm to compute the centre in stage 1; the value is either "OLS", corresponding to an ordinary linear regression under Gaussian assumptions to compute the mean, or "LAD", corresponding to a quantile regression to compute the median. The recommended default option is "LAD". For the quantile regression, the function calls \code{quantreg::rq} and the median is estimated using either the "br" (smaller samples) or "sfn" (larger samples and sparse problems) algorithm depending the sample size, for more details see \code{?quantreg::rq}.
#' @param transformed a logical indicating whether the quantitative response \code{Y} should be transformed using a rank-based method to resemble a normal distribution; recommended for traits with non-symmetric distribution. The default option is \code{FALSE}.
#' @param Xchr a logical indicator for whether the analysis is for X-chromosome SNPs.
#' @param origLev a logical indicator for whether the reported p-values should also include original Levene's test.
#' @param centre a character indicating whether the absolute deviation should be calculated with respect to "median" or "mean" in the traditional sex-specific and Fisher combined Levene's test p-values (three tests) for X-chromosome. The default value is "median". This option applies to sex-specific analysis using original Levene's test (i.e. when \code{regression}$$=$$\code{TRUE}).
#'
#' @param related optional: a logical indicating whether the samples should be treated as related; if \code{TRUE} while no relatedness covariance information is given, it is then estimated under a \code{cov.structure} and assumes this structure among all within-group errors pertaining to the same pair/cluster if specified using \code{clust}. This option currently only applies to autosomal SNPs.
#' @param cov.structure optional: should be one of standard classes of correlation structures listed in \code{corClasses} from \pkg{R} package \pkg{nlme}. See \code{?corClasses}. The most commonly used option is \code{corCompSymm} for a compound symmetric correlation structure. This option currently only applies to autosomal SNPs.
#' @param clust optional: a factor indicating the grouping of samples; it should have at least two distinct values. It could be the family ID (FID) for family studies. This option currently only applies to autosomal SNPs.
#'
#' @import stats
#' @import quantreg
#' @import methods
#' @import nlme
#'
#' @return a vector of Levene's test regression p-values according to the models
#' specified.
#' @export scaleReg
#' @note We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
#' the variance values tend to be proportional to mean values when stratified by \code{GENO}.
#'
#' @examples
#' N <- 1000
#' genoDAT <- rbinom(N, 2, 0.3)
#' sex <- rbinom(N, 1, 0.5)+1
#' Y <- rnorm(N)
#' covar <- matrix(rnorm(N*10), ncol=10)
#'
#' # vanilla example:
#' scaleReg(GENO=list(genoDAT, genoDAT), Y=Y, COVAR=covar)
#' scaleReg(GENO=list(genoDAT, genoDAT), Y=Y, COVAR=covar, genotypic=TRUE)
#' scaleReg(GENO=list(genoDAT, genoDAT), Y=Y, COVAR=covar, origLev = TRUE)
#' scaleReg(GENO=list(genoDAT, genoDAT), Y=Y, COVAR=covar, origLev = TRUE, SEX=sex)
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#' @references Deng WQ, Mao S, Kalnapenkis A, Esko T, Magi R, Pare G, Sun L. (2019) Analytical strategies to include the X-chromosome in variance heterogeneity analyses: Evidence for trait-specific polygenic variance structure. \emph{Genet Epidemiol}. \strong{43}(7):815-830. \doi{10.1002/gepi.22247}. PMID:31332826.
#' @references Gastwirth JL, Gel YR, Miao W. (2009). The Impact of Levene's Test of Equality of Variances on Statistical Theory and Practice." \emph{Statistical Science}. \strong{24}(3) 343-360, \doi{10.1214/09-STS301}.
#' @references Soave D, Sun L. (2017). A generalized Levene's scale test for variance heterogeneity in the presence of sample correlation and group uncertainty. \emph{Biometrics}. \strong{73}(3):960-971. \doi{10.1111/biom.12651}. PMID:28099998.




scaleReg <- function(GENO, Y, COVAR = NULL, SEX = NULL, Xchr = FALSE, transformed=FALSE, loc_alg = "LAD", related = FALSE, cov.structure = "corCompSymm", clust = NULL, genotypic = FALSE,  origLev = FALSE, centre = "median"){

  if (missing(Y))
      stop("The quantitative trait input is missing.")

  if (class(Y)!="numeric")
      stop("Please make sure the quantitaitve trait is a numeric vector.")

  if (missing(GENO))
     stop("The genotype input is missing.")

  if (!("matrix" %in% class(GENO) | "data.frame" %in% class(GENO) | "list" %in% class(GENO) | "vector" %in% class(GENO) | "integer" %in% class(GENO) | "numeric" %in% class(GENO))){
    stop("Please make sure the genotype data is an object of vector, matrix, data.frame for discrete genotypes or a list for dosage genotypes.")
  }


  if ("vector" %in% class(GENO) | "integer" %in% class(GENO) | "numeric" %in% class(GENO)){

  if (Xchr) {

      p_out <- tryCatch(leveneRegX_per_SNP(geno_one=GENO, SEX = SEX, Y=Y, COVAR = COVAR, transformed=transformed, loc_alg = loc_alg, genotypic=genotypic),  error=function(e) NA)

          if (origLev) {

          p_out2 <- tryCatch(leveneTests_per_SNP(geno_one = GENO, SEX = SEX, Y=Y, centre = "median", transformed=TRUE),  error=function(e) NA)

        if (is.null(names(GENO))){

        if (is.na(p_out2)){
          outputRes <- cbind(data.frame("CHR" = "X", "SNP" = "SNP", "gS" = p_out))
          } else {
          outputRes <- cbind(data.frame("CHR" = "X", "SNP" = "SNP", "gS" = p_out), p_out2)
          }
             } else {

        if (is.na(p_out2)){
          outputRes <-  cbind(data.frame("CHR" = "X", "SNP" = names(GENO), "gS" = p_out))
          } else {
        	outputRes <-  cbind(data.frame("CHR" = "X", "SNP" = names(GENO), "gS" = p_out), p_out2)
          }
             }

             } else {

        if (is.null(names(GENO))){

     	outputRes <- data.frame("CHR" = "X", "SNP" = "SNP", "gS" = p_out)

               } else {
    	outputRes <- data.frame("CHR" = "X", "SNP" = names(GENO), "gS" = p_out)
               }
          }

  } else {

      p_out <- tryCatch(leveneRegA_per_SNP(geno_one = GENO, Y = Y, COVAR = COVAR, transformed=transformed, loc_alg = loc_alg, related = related, cov.structure = cov.structure, clust = clust, genotypic = genotypic),  error=function(e) NA)

      if (origLev) {

        p_out2 <- tryCatch(leveneTests_per_SNP(geno_one = GENO, Y=Y, centre = "median", transformed=TRUE),  error=function(e) NA)

        if (is.null(colnames(GENO))|is.null(names(GENO))){
          outputRes <- cbind(data.frame("SNP" = "SNP", "gS" = p_out), p_out2)

        } else {
          outputRes <-  cbind(data.frame("SNP" = names(GENO), "gS" = p_out), p_out2)
        }

      } else {

        if (is.null(colnames(GENO))|is.null(names(GENO))){
          outputRes <- data.frame("SNP" = "SNP", "gS" = p_out)

        } else {

          outputRes <- data.frame("SNP" = colnames(GENO), "gS" = p_out)
        }
      }

      }

  }


  if ("data.frame"  %in% class(GENO) | "matrix" %in% class(GENO)){

    if (Xchr) {

      output_test <- apply(GENO, 2, function(ee) tryCatch(leveneRegX_per_SNP(geno_one=ee, SEX = SEX, Y=Y, COVAR = COVAR, transformed=transformed, loc_alg = loc_alg, genotypic=genotypic), error=function(e) NULL))

      p_out <- output_test
      change_to_NA <- sapply(output_test, is.null)

      if (sum(change_to_NA) == length(output_test)){
        p_out <- rep(NA, length(output_test))
      } else {
        NA_out <- rep(NA, length(p_out[[which(!change_to_NA)[1]]]))
        names(NA_out) <- names(p_out[[which(!change_to_NA)[1]]])
        p_out[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
      }

      if (is.null(dim(p_out))){
        p_out <- do.call(rbind, p_out)
      }

      if (origLev) {

        output_test2 <- apply(GENO, 2, function(ee) tryCatch(leveneTests_per_SNP(geno_one = ee, SEX = SEX, Y=Y, centre = centre, transformed=transformed), error=function(e) NULL))

        p_out2 <- output_test2
        change_to_NA <- sapply(output_test2, is.null)

        if (sum(change_to_NA) == length(output_test2)){
          p_out2 <- NA
        } else {
          NA_out <- rep(NA, length(p_out2[[which(!change_to_NA)[1]]]))
          names(NA_out) <- names(p_out2[[which(!change_to_NA)[1]]])
          p_out2[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
        }

        if (is.null(dim(p_out2)) & !is.na(p_out2)){
          p_out2 <- do.call(rbind, p_out2)
        }

        if (is.null(colnames(GENO))){

          if (sum(change_to_NA) == length(output_test2)) {
          outputRes <- cbind(data.frame("CHR" = "X", "SNP" = paste("SNP_",1:dim(GENO)[2],sep=""), "gS" = p_out))
          } else {
          outputRes <- cbind(data.frame("CHR" = "X", "SNP" = paste("SNP_",1:dim(GENO)[2],sep=""), "gS" = p_out), p_out2)
          }
        } else {
          if (sum(change_to_NA) == length(output_test2)) {
            outputRes <-  cbind(data.frame("CHR" = "X", "SNP" = colnames(GENO), "gS" = p_out))
          } else {
          outputRes <-  cbind(data.frame("CHR" = "X", "SNP" = colnames(GENO), "gS" = p_out), p_out2)
          }
        }

      } else {

        if (is.null(colnames(GENO))){

          outputRes <- data.frame("CHR" = "X", "SNP" = paste("SNP_",1:dim(GENO)[2],sep=""), "gS" = p_out)

        } else {

          outputRes <- data.frame("CHR" = "X", "SNP" = colnames(GENO), "gS" = p_out)
        }
      }

    } else {

      output_test <- apply(GENO, 2, function(ee) tryCatch(leveneRegA_per_SNP(geno_one = ee, Y = Y, COVAR = COVAR, transformed=transformed, loc_alg = loc_alg, related = related, cov.structure = cov.structure, clust = clust, genotypic = genotypic), error=function(e) NULL))

      p_out <- output_test
      change_to_NA <- sapply(output_test, is.null)

      if (sum(change_to_NA) == length(output_test)){
        p_out <- rep(NA, length(output_test))
      } else {
        NA_out <- rep(NA, length(p_out[[which(!change_to_NA)[1]]]))
        names(NA_out) <- names(p_out[[which(!change_to_NA)[1]]])
        p_out[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
      }

      if (is.null(dim(p_out))){
        p_out <- do.call(rbind, p_out)
      }

      if (origLev) {

        output_test2 <- apply(GENO, 2, function(ee) tryCatch(leveneTests_per_SNP(geno_one = ee, SEX = SEX, Y=Y, centre = centre, transformed=transformed), error=function(e) NULL))

        p_out2 <- output_test2
        change_to_NA <- sapply(output_test2, is.null)

        if (sum(change_to_NA) == length(output_test2)){
          p_out2 <- NA
        } else {
          NA_out <- rep(NA, length(p_out2[[which(!change_to_NA)[1]]]))
          names(NA_out) <- names(p_out2[[which(!change_to_NA)[1]]])
          p_out2[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
        }

        if (is.null(dim(p_out2)) & !is.na(p_out2)){
          p_out2 <- do.call(rbind, p_out2)
        }

        if (is.null(colnames(GENO))){

          if (sum(change_to_NA) == length(output_test2)){
            outputRes <- cbind(data.frame("SNP" = paste("SNP_",1:dim(GENO)[2],sep=""), "gS" = p_out))
          } else {
          outputRes <- cbind(data.frame("SNP" = paste("SNP_",1:dim(GENO)[2],sep=""), "gS" = p_out), p_out2)
          }

        } else {
          if (sum(change_to_NA) == length(output_test2)){
          outputRes <-  cbind(data.frame("SNP" = colnames(GENO), "gS" = p_out))
          } else {
          outputRes <-  cbind(data.frame("SNP" = colnames(GENO), "gS" = p_out), p_out2)
          }

        }

      } else {

        if (is.null(colnames(GENO))){
          outputRes <- data.frame("SNP" = paste("SNP_",1:dim(GENO)[2],sep=""), "gS" = p_out)

        } else {
          outputRes <- data.frame("SNP" = colnames(GENO), "gS" = p_out)
        }
      }

    }
  }


  if ("list" %in% class(GENO)){

      if (Xchr) {

        output_test <- lapply(GENO, function(ee) tryCatch(leveneRegX_per_SNP(geno_one=ee, SEX = SEX, Y=Y, COVAR = COVAR, transformed=transformed, loc_alg = loc_alg, genotypic=genotypic), error=function(e) NULL))

        p_out <- output_test
        change_to_NA <- unlist(lapply(output_test, is.null))

        if (sum(change_to_NA) == length(output_test)){
          p_out <- list(NA)[rep(1, length(output_test))]
        } else {
          NA_out <- rep(NA, length(p_out[[which(!change_to_NA)[1]]]))
          names(NA_out) <- names(p_out[[which(!change_to_NA)[1]]])
          p_out[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
        }

      if (origLev) {

      output_test2 <- lapply(GENO, function(ee) tryCatch(leveneTests_per_SNP(geno_one = ee, SEX = SEX, Y=Y, centre = centre, transformed=transformed), error=function(e) NULL))

      p_out2 <- output_test2
      change_to_NA <- unlist(lapply(output_test2, is.null))

      if (sum(change_to_NA) == length(output_test2)){
        p_out2 <- NA
      } else {
        NA_out <- rep(NA, length(p_out2[[which(!change_to_NA)[1]]]))
        names(NA_out) <- names(p_out2[[which(!change_to_NA)[1]]])
        p_out2[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
      }

        if (is.null(names(GENO))){

          if (sum(change_to_NA) == length(output_test2)){
            outputRes <- cbind(data.frame("CHR" = "X", "SNP" = paste("SNP_",1:length(GENO), sep=""), "gS" = unlist(p_out)))
          } else {
          outputRes <- cbind(data.frame("CHR" = "X", "SNP" = paste("SNP_",1:length(GENO), sep=""), "gS" = unlist(p_out)), do.call(rbind, p_out2))
          }

        } else {

          if (sum(change_to_NA) == length(output_test2)){
            outputRes <-  cbind(data.frame("CHR" = "X", "SNP" = names(GENO), "gS" = unlist(p_out)))
          } else {
          outputRes <-  cbind(data.frame("CHR" = "X", "SNP" = names(GENO), "gS" = unlist(p_out)), do.call(rbind, p_out2))
          }
        }

      } else {

        if (is.null(names(GENO))){
          outputRes <- data.frame("CHR" = "X", "SNP" = paste("SNP_",1:length(GENO), sep=""), "gS" = unlist(p_out))

        } else {
          outputRes <- data.frame("CHR" = "X", "SNP" = names(GENO), "gS" = unlist(p_out))
        }
      }

    } else {

      output_test <- lapply(GENO, function(ee) tryCatch(leveneRegA_per_SNP(geno_one = ee, Y = Y, COVAR = COVAR, transformed=transformed, loc_alg = loc_alg, related = related, cov.structure = cov.structure, clust = clust, genotypic = genotypic), error=function(e) NULL))

      p_out <- output_test
      change_to_NA <- unlist(lapply(output_test, is.null))

      if (sum(change_to_NA) == length(output_test)){
        p_out <- list(NA)[rep(1, length(output_test))]
      } else {
        NA_out <- rep(NA, length(p_out[[which(!change_to_NA)[1]]]))
        names(NA_out) <- names(p_out[[which(!change_to_NA)[1]]])
        p_out[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
      }

      if (origLev) {

        output_test2 <- lapply(GENO, function(ee) tryCatch(leveneTests_per_SNP(geno_one = ee, SEX = SEX, Y=Y, centre = centre, transformed=transformed), error=function(e) NULL))

        p_out2 <- output_test2
        change_to_NA <- unlist(lapply(output_test2, is.null))

        if (sum(change_to_NA) == length(output_test2)){
          p_out2 <- NA
        } else {
          NA_out <- rep(NA, length(p_out2[[which(!change_to_NA)[1]]]))
          names(NA_out) <- names(p_out2[[which(!change_to_NA)[1]]])
          p_out2[which(change_to_NA)] <- list(NA_out)[rep(1, sum(change_to_NA))]
        }

        if (is.null(names(GENO))){

          if (sum(change_to_NA) == length(output_test2)){
            outputRes <- cbind(data.frame("SNP" = paste("SNP_",1:length(GENO), sep=""), "gS" = unlist(p_out)))
          } else {
            outputRes <- cbind(data.frame("SNP" = paste("SNP_",1:length(GENO), sep=""), "gS" = unlist(p_out)), do.call(rbind, p_out2))
          }

        } else {
          if (sum(change_to_NA) == length(output_test2)){
            outputRes <-  cbind(data.frame("SNP" = names(GENO), "gS" = unlist(p_out)))
          } else {
          outputRes <-  cbind(data.frame("SNP" = names(GENO), "gS" = unlist(p_out)), do.call(rbind, p_out2))
          }
        }

      } else {

        if (is.null(names(GENO))){
          outputRes <- data.frame("SNP" = paste("SNP_",1:length(GENO), sep=""), "gS" = unlist(p_out))

        } else {
          outputRes <- data.frame("SNP" = names(GENO), "gS" = unlist(p_out))
        }
      }

    }

}

  rownames(outputRes) <- NULL

  return(outputRes)
 }

