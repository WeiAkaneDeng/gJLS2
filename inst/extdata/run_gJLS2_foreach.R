#!/usr/bin/env Rscript

if("optparse" %in% rownames(installed.packages()) == FALSE) {
print("optparse not installed, trying to intall now ...")
install.packages("optparse", repos='http://cran.us.r-project.org')
}

require("optparse")
 
option_list = list(
  make_option(c("-b", "--bfile"), type="character", default=NULL, 
              help="genotype dataset file name", metavar="character"),
  make_option(c("-p", "--pfile"), type="character", default=NULL, 
              help="pheno and covariate dataset file name", metavar="character"),
  make_option(c("-m", "--pheno"), type="character", default=NULL, 
              help="phenotype name in pheno file", metavar="vector"),
  make_option(c("-v", "--covar"), type="character", default=NULL, 
              help="covariate name in pheno file", metavar="vector"),
  make_option(c("-q", "--center"), type="character", default="median", 
              help="center option in gJLS2", metavar="character"),
  make_option(c("-g", "--genotypic"), type="logical", default="TRUE", 
              help="genotypic option in gJLS2 for scale analysis", metavar="character"),
  make_option(c("-t", "--transform"), type="logical", default="FALSE", 
              help="phenotype transform option in gJLS2", metavar="character"),
  make_option(c("-x", "--Xchr"), type="logical", default="FALSE", 
              help="Xchr option in gJLS2", metavar="character"),
  make_option(c("-n", "--nThreads"), type="integer", default="1", 
              help="number of Threads used", metavar="integer"),
  make_option(c("-a", "--analysis"), type="integer", default="2", 
              help="location (0), scale (1), or joint-location-scale analysis (2)", 
              metavar="integer"),
  make_option(c("-w", "--write"), type="integer", default= 50, 
              help="chunk size to write in output file; /n
                    default is 50 lines", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-s", "--sumfile"), type="character", default=NULL, 
              help="summary statistics dataset file name;/n 
                    must contain column names SNP, gL and gS", metavar="character")
); 
 

opt_parser <-OptionParser(option_list=option_list)
arguments <- parse_args(opt_parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args


nThread <- opt$nThreads
maxCores <- parallel::detectCores(logical = TRUE)
cl <- parallel::makeCluster(min(maxCores-1, nThread))
doParallel::registerDoParallel(cl)


if("parallel" %in% rownames(installed.packages()) == FALSE) {
cat("parallel not installed, trying to intall now ...")
install.packages("parallel", repos='http://cran.us.r-project.org')
}
require("parallel")	


if("doParallel" %in% rownames(installed.packages()) == FALSE) {
cat("doParallel not installed, trying to intall now ...")
install.packages("doParallel", repos='http://cran.us.r-project.org')
}
require("doParallel")	


if("iterators" %in% rownames(installed.packages()) == FALSE) {
cat("iterators not installed, trying to intall now ...")
install.packages("iterators", repos='http://cran.us.r-project.org')
}
require("iterators")	



if("foreach" %in% rownames(installed.packages()) == FALSE) {
cat("foreach not installed, trying to intall now ...")
install.packages("foreach", repos='http://cran.us.r-project.org')
}
library("foreach")	


if (is.null(opt$sumfile)){

bimf <- opt$bfile
phenof <-opt$pfile

phenoNames <- opt$pheno
covarNames <- strsplit(opt$covar, ",")[[1]]
chunk_size <- opt$write
cat(paste("Writing in chunk size of", chunk_size, "\n"))

## additional options:

centre <- opt$center; 
cat(paste("Using center option", centre, "\n"))
genotypic <- opt$genotypic
cat(paste("Using genotypic option", genotypic, "\n"))
transformY <- opt$transform
cat(paste("Using transform option", transformY, "\n"))
xchr <- opt$Xchr
cat(paste("Using Xchr option", xchr, "\n"))
out <- opt$out
runA <- opt$analysis
cat(paste("Running", c("location only", "scale only", "joint-location and scale")[runA+1], "analysis", "\n"))


if("gJLS2" %in% rownames(installed.packages()) == FALSE) {
cat("gJLS2 not installed, trying to intall now ...")
#install.packages("gJLS2", repos='http://cran.us.r-project.org')
library("devtools")
devtools::install_github("WeiAkaneDeng/gJLS2")
}

require("gJLS2")

## checking pheno file

if("BEDMatrix" %in% rownames(installed.packages()) == FALSE) {
print("BEDMatrix not installed, trying to intall now ...")
install.packages("BEDMatrix", repos='http://cran.us.r-project.org')
}

if("BGData" %in% rownames(installed.packages()) == FALSE) {
print("BGData not installed, trying to intall now ...")
install.packages("BGData", repos='http://cran.us.r-project.org', dependencies=T)
}

## checking inputs to be bed, fam, bim files

require("BGData")
require("BEDMatrix")

bedFiles <- BEDMatrix(bimf)

cat(paste("linking phenotype file", phenof, "\n"))

bg <- as.BGData(bedFiles, alternatePhenotypeFile = paste0(phenof))
	
## CHECKING ALL INPUT FILES AGAIN:

pheno_dat <- pheno(bg)
geno_dat <- geno(bg)

if (sum(grepl("SEX", names(pheno_dat)))>1){
	names(pheno_dat)[grepl("SEX", names(pheno_dat))] <- c("SEX", paste("SEX", 1:(dim(pheno_dat[grepl("SEX", names(pheno_dat))])[2]-1), sep=""));
}


if (!is.null(covarNames)){

if (sum(grepl("sex|SEX|Sex", covarNames)) > 0){

	SEX_cov <- as.integer(pheno_dat[,names(pheno_dat) %in% covarNames][grepl("sex|SEX|Sex", covarNames)][,1])
	covarNames_use <- covarNames[!grepl("sex|SEX|Sex", covarNames)]
	SEX_cov_PLINK <- ifelse(SEX_cov==0, 2, SEX_cov)

cat(paste("Covariates include", covarNames, " from covariate/pheno file \n"))

} else{

	SEX_cov <- pheno_dat$SEX
	SEX_cov_PLINK <- ifelse(SEX_cov==0, 2, SEX_cov)
	covarNames_use <- covarNames

cat(paste("Covariates did not include SEX, taking SEX from .fam file\n"))
}

cat(paste("Writing results to output", out, "\n"))

## writing results by chunks of 50 SNPs to avoid loss in interruption

iteraR <- round(dim(geno_dat)[2]/chunk_size)

if (runA == 0) {

chunked_df <- iterators::iter(geno_dat[,(1):min(dim(geno_dat)[2], chunk_size)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages= "foreach") %dopar% {
gJLS2::locReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY)
}

	
} else if (runA == 1) {

chunked_df <- iterators::iter(geno_dat[,(1):min(dim(geno_dat)[2], chunk_size)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::scaleReg(GENO = to_compute,Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY, genotypic = genotypic, centre = centre)
}

	
} else {

chunked_df <- iterators::iter(geno_dat[,(1):min(dim(geno_dat)[2], chunk_size)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::gJLS2(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr, transformed = transformY, genotypic = genotypic, centre = centre)
}

}

write.table(final_output, file = out, col.names=T, row.names=F, quote=F, sep="\t")


if (iteraR > 1) {
	
for (j in 2:iteraR){
	
	cat(paste("Running the", j, "th chunk", "\n"))

if (j == iteraR){
	
if (runA == 0) {	

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(iteraR-1)):dim(geno_dat)[2]], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::locReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY)
}

} else if (runA == 1) {


chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(iteraR-1)):dim(geno_dat)[2]], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::scaleReg(GENO = to_compute,Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY, genotypic = genotypic, centre = centre)
}

} else {

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(iteraR-1)):dim(geno_dat)[2]], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::gJLS2(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY, genotypic = genotypic, centre = centre)
}
}

write.table(final_output, file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")

} else {	

if (runA == 0) {	

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::locReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY)
}

} else if (runA == 1) {

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::scaleReg(GENO = to_compute,Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY, genotypic = genotypic, centre = centre)
}


} else {

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::gJLS2(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr,  transformed = transformY, genotypic = genotypic, centre = centre)
}

}
write.table(final_output, file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
}
}
}



} else {


cat(paste("Writing results to output", out, "\n"))
	
iteraR <- round(dim(geno_dat)[2]/chunk_size)

if (runA == 0) {	

chunked_df <- iterators::iter(geno_dat[,(1):(chunk_size)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::locReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY)
}

} else if (runA == 1) {

chunked_df <- iterators::iter(geno_dat[,(1):(chunk_size)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::scaleReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY, genotypic = genotypic, centre = centre)
}


} else {

chunked_df <- iterators::iter(geno_dat[,(1):(chunk_size)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::gJLS2(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY, genotypic = genotypic, centre = centre)
}

}

write.table(final_output, file = out, col.names=T, row.names=F, quote=F, sep="\t")

if (iteraR > 1) {

for (j in 2:iteraR){

	cat(paste("Running the", j, "th chunk", "\n"))

if (j == iteraR){
	
if (runA == 0) {	

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(iteraR-1)):dim(geno_dat)[2]], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::locReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY)
}

} else if (runA == 1) {

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(iteraR-1)):dim(geno_dat)[2]], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::scaleReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY, genotypic = genotypic, centre = centre)
}


} else {

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(iteraR-1)):dim(geno_dat)[2]], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::gJLS2(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY, genotypic = genotypic, centre = centre)
}


}

write.table(final_output, file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")


	} else {	

if (runA == 0) {	

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::locReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY)
}


} else if (runA == 1) {

chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::scaleReg(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY, genotypic = genotypic, centre = centre)
}

} else {


chunked_df <- iterators::iter(geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], by="column", chunksize=round(sqrt(chunk_size))) 

final_output <- foreach(to_compute=chunked_df, .combine = rbind, .packages="foreach") %dopar% {
gJLS2::gJLS2(GENO = to_compute, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, nCores= nThread, transformed = transformY, genotypic = genotypic, centre = centre)
}




}

write.table(final_output, file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")


	}	
	}
}

}


stopCluster(cl)

} else {

	out <- opt$out
	sum_file <- read.table(opt$sumfile, head=T)	
	final_output <- gJLS2::gJLS2s(gL = sum_file$gL, gS = sum_file$gS)
	sum_file_out <- sum_file
	sum_file_out$gJLS <- formatC(final_output$gJLS, digits=3, format="e")
	write.table(sum_file_out, file = out, col.names=T, row.names=F, quote=F, sep="\t")	
}






		