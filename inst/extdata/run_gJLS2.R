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
  make_option(c("-n", "--nTasks"), type="integer", default="1", 
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

library("gJLS2")

if("BEDMatrix" %in% rownames(installed.packages()) == FALSE) {
print("BEDMatrix not installed, trying to intall now ...")
install.packages("BEDMatrix", repos='http://cran.us.r-project.org')
}

if("BGData" %in% rownames(installed.packages()) == FALSE) {
print("BGData not installed, trying to intall now ...")
install.packages("BGData", repos='http://cran.us.r-project.org', dependencies=T)
}


if("bigstatsr" %in% rownames(installed.packages()) == FALSE) {
print("bigstatsr not installed, trying to intall now ...")
install.packages("bigstatsr", repos='http://cran.us.r-project.org')
}


nTasks <- opt$nTasks
nMaxcores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
nTasks_use <- as.numeric(max(1, min(nTasks, nMaxcores, na.rm=T), na.rm=T))
cat(paste("number of cores initiated is", nTasks_use, "\n"))

if (nTasks_use > 1){
if("parallel" %in% rownames(installed.packages()) == FALSE) {
cat("parallel not installed, trying to intall now ...")
install.packages("parallel", repos='http://cran.us.r-project.org')
}
library(parallel)	
}

checkForTcltk <- function(){
    if ("tcltk" %in% loadedNamespaces()){
        warning("This function cannot be used because the R tcltk package is loaded. Changing to the default number of cores.")
    nTasks_use = 1
    }
}
checkForTcltk()


### write generic functions for each analysis:

if (runA == 0) {

runFunction <- function(ee){
gJLS2::locReg(GENO = ee, Y = Y_PLINK, SEX = SEX_cov_PLINK, COVAR = COV_plink, Xchr=xchr, transformed = transformY)
}

} else if (runA == 1) {

runFunction <- function(ee){
gJLS2::scaleReg(GENO = ee, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = COV_plink, Xchr=xchr, transformed = transformY, genotypic = genotypic, centre = centre)
}

} else {
runFunction <- function(ee){
gJLS2::gJLS2(GENO = ee, Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = COV_plink, Xchr=xchr, transformed = transformY, genotypic = genotypic, centre = centre)
}
}


#### checking inputs to be bed, fam, bim files

require("BGData")
require("BEDMatrix")
bedFiles <- BEDMatrix(bimf)


## writing results by chunks of 100 SNPs to avoid loss in interruption
## make sure chunk size is even
nbSNPs <- dim(bedFiles)[2]

cat(paste("Reading", nbSNPs, "SNPs from the bed file \n"))

chunk_size <- min(chunk_size, nbSNPs)

iteraR <- max(1, ceiling(nbSNPs/chunk_size))

cat(paste("Expecting to run", iteraR, "chunks in total \n"))


cat(paste("linking phenotype file", phenof, "\n"))
bg <- as.BGData(bedFiles, alternatePhenotypeFile = paste0(phenof))
pheno_dat <- pheno(bg)

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

COV_plink <- pheno_dat[,names(pheno_dat) %in% covarNames_use]

} else {

COV_plink <- NULL
SEX_cov_PLINK <- as.numeric(as.character(read.table(bimf)[,5]))
SEX_cov_PLINK[SEX_cov_PLINK==-9] <- NA

cat(paste("Trying to include SEX from .fam file \n"))
}
	

if (iteraR == 1) {

geno_dat <- geno(bg)

cat(paste("Running the 1st chunk", "\n"))

final_output <- chunkedMap(X = geno_dat, FUN = runFunction, nCores=nTasks_use)

cat(paste("Writing results to output", out, "\n"))

write.table(final_output, file = out, col.names=T, row.names=F, quote=F, sep="\t")


} else {
		
chunk_list <- bigstatsr:::CutBySize(nbSNPs, nb = chunk_size)

for (j in 1:iteraR){
	
	cat(paste("Running chunk num", j, "\n"))
	geno_dat <- geno(bg)[,chunk_list[j,]]
	final_output <- chunkedMap(X = geno_dat, FUN = runFunction, nCores=nTasks_use)
	cat(paste("Writing results to output", out, "\n"))

if (j==1){
	write.table(final_output, file = out, col.names=T, row.names=F, quote=F, sep="\t")
} else {
	write.table(final_output, file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
}

}
}


} else {
	out <- opt$out
	sum_file <- read.table(opt$sumfile, head=T)	
	final_output <- gJLS2::gJLS2s(gL = sum_file$gL, gS = sum_file$gS)
	sum_file_out <- sum_file
	sum_file_out$gJLS <- formatC(final_output$gJLS, digits=3, format="e")
	write.table(sum_file_out, file = out, col.names=T, row.names=F, quote=F, sep="\t")	
}






		