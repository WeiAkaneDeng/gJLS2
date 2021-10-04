#' generalized Joint-Location-Scale (gJLS) test with summary statistics
#'
#' This function takes as input the gL and gS p-values for each SNP and combine to produce the gJLS p-values. It is used for genome-wide analysis where only the gL or gS p-values are available, caution should be exercised when combing gL and gS p-values obtained from separate datasets.
#'
#' @param gL a vector of location p-values or a data.frame containing column names "SNP" and "gL".
#' @param gS a vector of scale p-values or a data.frame containing column names "SNP" and "gS".
#'
#' @importFrom methods is
#' @importFrom stats complete.cases
#' @importFrom stats na.exclude
#'
#' @return a vector of combined gJLS p-values for each SNP.
#' @export gJLS2s
#'
#' @note For a genome-scan, we recommend to run this in PLINK via the plugin function \code{gJLSPLINK}, especially for large datasets and those with more than 20 covariates.
#' @note We highly recommend to quantile-normally transform \code{Y} for non-symmetrically distributed traits. This is typically done to avoid ‘scale-effect’ when the variance values tend to be proportional to mean values when stratified by \code{GENO}, as observed by Pare et al. (2010) and Yang et al. (2011).
#' @note For the moment, only quantitative trait \code{Y} is accepted as the subsequent generalized joint location scale (gJLS) analyses require the variance be calculated on quantitative traits. However, we are working on to include binary response for the generalized JLS analyses in the next update of gJLS.
#'
#' @examples
#' gL <- data.frame("SNP" = paste("rs", 1:100, sep=""), "gL"=runif(100))
#' gS <- runif(100)
#'
#' gJLS2s(gL = gL, gS=gS)
#'
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}, Lei Sun \email{sun@utstat.toronto.edu}
#'
#' @references Soave D, Corvol H, Panjwani N, Gong J, Li W, Boëlle PY, Durie PR, Paterson AD, Rommens JM, Strug LJ, Sun L. (2015). A Joint Location-Scale Test Improves Power to Detect Associated SNPs, Gene Sets, and Pathways. \emph{American Journal of Human Genetics}. 2015 Jul 2;\strong{97}(1):125-38. \doi{10.1016/j.ajhg.2015.05.015}. PMID: 26140448; PMCID: PMC4572492.



gJLS2s <- function(gL, gS){

  if (is.null(dim(gL)) & is.null(dim(gS))){

    if (length(gL)!= length(gS))
      stop("Please make sure the two input p-value vectors are of the same length.")

    merged_dat <- data.frame("gL" = gL, "gS" = gS)
    merged_dat$gJLS <- 1-pchisq(-2*log(merged_dat$gL)-2*log(merged_dat$gS), 4)

  } else {

    if (!is.null(dim(gL)) & !is.null(dim(gS))) {

      gL <- as.data.frame(gL)
      gS <- as.data.frame(gS)

      if (dim(gL)[1]!=dim(gS)[1])
        stop("Please make sure the two input p-value data.frames are of the same length.")

      if ("SNP" %in% names(gL) & "gL" %in% names(gL))
        stop("gL is a data.frame but does not contain the column names 'SNP' and 'gL'.")

      if ("SNP" %in% names(gS) & "gS" %in% names(gS))
        stop("gS is a data.frame but does not contain the column names 'SNP' and 'gL'.")

      suppressMessages(merged_dat <- plyr::join(gL, gS))
      merged_dat$gJLS <- 1-pchisq(-2*log(merged_dat$gL)-2*log(merged_dat$gS), 4)

    } else if (!is.null(dim(gL))){

      if (dim(gL)[1]!= length(gS))
        stop("Please make sure the two input p-value data.frame/vectors are of the same length.")

      merged_dat <- cbind(gL, "gS" = gS)
      merged_dat$gJLS <- 1-pchisq(-2*log(merged_dat$gL)-2*log(merged_dat$gS), 4)

    } else {

      if (length(gL)!= dim(gS)[1])
        stop("Please make sure the two input p-value data.frame/vectors are of the same length.")
      merged_dat <- cbind("gL" = gL, gS)
      merged_dat$gJLS <- 1-pchisq(-2*log(merged_dat$gL)-2*log(merged_dat$gS), 4)

    }
}

  rownames(merged_dat) <- NULL

  return(merged_dat)
}


