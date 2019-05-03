 
## Repository for the "Genome-wide nucleotide patterns and potential mechanisms for genome divergence following domestication in maize and soybean" project. ##

**Note: in this depository we did not store the original data mainly because there are already available online in the original publication. Here we will document the analysis in detail and present a clear pipeline.**

### Outline of the repository ###


In order to better guide the visitors about the repository, here we briefly introduce the outline of the repository
 
**1. Filtering common SNP set** 

-Website address for the original data set

-Script for filtering common SNP set from the original dataset


**2. Genome-wide base composition** 

-Script for base composition calculation from genome-wide SNPs

-Dataset of base composition value for each accession

**3. Base composition among substitution types** 

-Script for calculation of base-composition conditional on each substitution types

-Dataset of base composition value conditional on each substitution type 

**4. Base composition distribution at different genomic regions**
  
-Predict SNP effect with SnpEff

-Script to classify SNPs in to different genomic annotation sets based on their predicted SNP effect 

-Script to calculate the base composition from SNPs in different genomic regions

-Script to randomly sample the same amount of snps from both genic and nongenic to check the snp_effect

-Script to calculate the base composition distribution across segments along chromsomes 

-Script to calculate the minor allele frequency distribution

-Script to calculate base composition at nongenic and genic conditional on pericentromeric, and chromosome separately

-Script to calculate the base composition distribution from selection sweep and non-selective sweep regions

-Data for percentage of SNPs in each of the genomic annotation sets

-Data for base composition distribution for different genomic annotation sets

-Data for base composition along chromosome (only one chromosome)

-Data for base composition for genic nongenic conditional on peri nonperi

-Data for base composition at selective sweep and nonsweep region

**5. Motif enrichment analysis**

-Script to calculate the frequency of 96 motifs from common SNP set at SNP site

-Script to calculate the frequency of 96 motifs from commmon SNP set at random site

-Script to calculate 4 solar-UV related motif genic, nongenic, peri, nonperi, conditional on methylated and unmethylated regions

-Data for motifs from random site and SNP site used for the motif enrichment anlaysis

-Data for motifs at different genomic regions

**6. Mutation spectrum from population-private SNPs**

-Blast to *Medicago truncatula* for the orthlogous regions to identify the ancestral state of the soybean SNPs

-Scripts for filtering population private SNPs from the original dataset

-Script for caculating frequency of population private SNPs

-Data for conducting mutation-spectrum analysis 

-Data for pop-private motif across SNP bins

**7. GWAS for base-composition across polymorphic sites**

-UV damage repair gene list in maize and soybean

-Data for regional Manhattan plot

-Data for haplotype analayis

