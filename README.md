 
## Repository for the "Genome-wide nucleotide patterns and potential mechanisms for genome divergence following domestication" project. ##

**Note: in this depository we did not store the original data mainly because there are already available online in the original publication. Here we will documenting the analysis in detail and make the pipeline clear.**

### Outline of the repository ###


In order to better guide the visitors about the repository, here we briefly introduce the outline of the repository
 
**1. Filtering common SNP set** 

-Website address for the original data set

-Script for filtering common SNP set from the original dataset


**2. Genome-wide base composition** 

-Script for calculation from genome-wide SNPs

-Dataset of base composition value for each accession

**3. Base composition among substitution types** 

-Script for calculation of base-composition conditional on each substitution types

-Dataset of base composition value conditional on each substitution type 

-R code for plotting Figure2

**4. Base composition distribution at different genomic regions**
  
-Predict SNP effect with SnpEff

-Classify SNPs in to different genomic annotation sets based on their predicted SNP effect 

-Script to calculate the base composition from SNPs in genic and nongenic regions

-Script to calculate the base composition distribution along chromsomes 

-Data to plot the Figure3

-R code to plot the Figure3

**5. Motif enrichment analysis**

-Script to calculate the frequency of 96 motifs from common SNP set at SNP site

-Script to calculate the frequency of 96 motifs from commmon SNP set at random site

-Script to calculate the motif frequency at genic and nongenic region

-Script to calculate the motif frequency at pericentromeric and nonpericentromeric region

-Script to differentiate SNPs into methylated and unmethylated regions

-Script to calculate the motif frequency from genic and nongenic SNPs conditional on methylated and unmethylated regions

-Data for conducting enrichment test and plotting Figure4 

**6. Mutation spectrum from population-private SNPs**

-Script for filtering population private SNPs from the original dataset

-Data for conducting mutation-spectrum analysis and plotting

-R code to plot Figure5

**7. GWAS for base-composition across polymorphic sites**

-UV damage repair gene list in maize and soybean

-Gene enrichment test

-Data for regional Manhattan plot

-Data for haplotype analayis

-R code including haplotype analysis and plotting Figure6