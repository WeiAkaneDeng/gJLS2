
# User (input data) scenarios

The gJSL2 analyzes each SNP at a time, it is straightforward to divide-and-conquer irrespective of your GUI choice. Here I list some possible user scenarios and our recommendation on which approach to use. We also report some performance metric to help you gauge the most appropriate option for your own computing specs.


- Assume the genotype data are either binary: **plink.bed**, **plink.bim**, and **plink.fam** or raw genotypes **plink.raw**; 

- Assume the phenotype data **pheno.txt** has 2+k columns where k is the number of phenotypes (including covariates if necessary) and the first two columns are "FID" and "IID" (which can be used to link to genotype files); the file should contain headers;

- Assume the summary file is **gJL2_summary_stats.txt** with columns **CHR** (not essential if SNP is present), **SNP**, **BP** (not essential if SNP is present), **gL**, **gS**.

Following each of the user scenarios, we provide some sample codes to get you started:


## **gJLS** with **INDVIDUAL-LEVEL** genotype and phenotype data.

  + *I have a relatively small dataset with < 5,000 samples and want to investigate particular markers or gene-set analysis (< 100).* <span style="color:green">The R GUI should be sufficient for this purpose. PLINK binary data can be read with "BEDMatrix" and the .raw genotypes can be read in directly as a data.frame.</span>
  
  
  + *I have a relatively small dataset with < 5,000 samples and want to perform a genome-wide analysis.* <span style="color:green">We recommend either the Rscript or PLINK R plugin as the analysis is straightforward and simple to break the job by chromosome or apply other filters within PLINK.</span>
  
  
  + *I have datasets with > 5,000 samples.* <span style="color:green">The Rscript option will work nicely on binary files using "BEDMatrix", the multiple cores option can speed up the computational time, and the user specified write size makes sure that no results are lost in the process.</span>


## **gS** with **INDVIDUAL-LEVEL** genotype and phenotype data.

- *I have "**INDVIDUAL-LEVEL**" genotype and phenotype data, but already ran the location analysis (GWAS) using other methods (e.g. PLINK or BOLT-LMM) or plan to use published GWAS p-values. I am interested in the "**gS**" analysis and the combined "**gJLS**" analyses using existing location p-values in place of the gL.* <span style="color:green">The user can run scale analysis alone using the preferred approach depending on sample size and computational requirements. See next scenario for combining gL and gS. </span>
  
  
## **gJLS** with **SUMMARY-LEVEL** genotype and phenotype data.

*I only have "**SUMMARY-LEVEL**" data (e.g. p-values) for location and scale association from external sources (e.g. GIANT) and am interested in the "**gJLS**" analysis.* <span style="color:green">Once the gL and gS *p*-values have been obtained, the function *gJLS2s()* can be called in R GUI via the Rscript to produce the combined gJLS *p*-values without access to individual-level data.</span>
  
  
