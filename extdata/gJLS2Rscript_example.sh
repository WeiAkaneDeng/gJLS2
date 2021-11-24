## R scripts for genotype/pheno data:


Rscript run_gJLS2.R --bfile ./input/chrX_5_snp \
--pfile ./input/Pheno.txt \
--pheno pheno1 \
--Xchr TRUE \
--nThreads 2 \
--covar SEX,covar1,covar2,covar3 \
--out ./output/testRun.results.txt





## R scripts for summary stats:

Rscript run_gJLS2.R --sumfile ./input/GIANT_BMI_chr16_gJLS_summary.txt \
--out ./output/GIANT_BMI_Sum.chr16_results.txt &

Rscript run_gJLS2.R --sumfile ./input/GIANT_Height_chr16_gJLS_summary.txt \
--out ./output/GIANT_Height_Sum.chr16_results.txt &