## This file contains the steps carried out for base composition distribution at different genomic regions##

**(Note: Maize part and soybean part are very similar, here mainly use maize part as an example to show the process)**

**1. Predict snp effect with snpEff**

- convert the hmp file to vcf file which will be used for the snpEff input
```
$ nohup perl Convert_from_hmp_to_vcf4_snpEff.pl 1 10 &
```

- Use snpEff to annotate the SNPs

1) Install the latest version of snpEff

```
$ wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
```

```
$ unzip snpEff_latest_core.zip
```

2) Build the database for MaizeZmB73_AGPV4 as the version4 data base are not pre-built in snpEff.

**(This part maize and soybean are different)**

**For maize:**

Download the genome annotation file: 

```
$ wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/gff3/zea_mays/Zea_mays.AGPv4.39.chr.gff3.gz
```

```
$ gunzip Zea_mays.AGPv4.39.chr.gff3.gz
```

```
$ mv Zea_mays.AGPv4.39.chr.gff3 genes.gff
```


Download the fasta file: 

```
$ wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
```

```
$ mv ZmB73_RefGen_V4.fa maizeZmB73_AGPV4.fa
```

Create database

```
$ java -jar snpEff.jar build -gff3 -v maizeZmB73_AGPV4 
```

**For soybean:**

Download the genome annotation file: 

```
$ wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/gff3/glycine_max/Glycine_max.Glycine_max_v2.0.39.chr.gff3.gz
```

```
$ gunzip Glycine_max.Glycine_max_v2.0.39.chr.gff3.gz
```

```
$ mv Glycine_max.Glycine_max_v2.0.39.chr.gff3 genes.gff
```


Download the fasta file: 

```
$ wget ftp://ftp.ensemblgenomes.org/pub/plants/release-39/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.0.dna.toplevel.fa.gz
```

```
$ mv Glycine_max.Glycine_max_v2.0.dna.toplevel.fa gmaxW82_v2.fa
```
Create database

```
$java -jar snpEff.jar build -gff3 -v gmaxW82_v2 
```

3) Run snpEff to annotate the converted vcf file

```
$ cd /snpEff/
```

**For maize:**

```
$ for i in {1..10}; do java -Xmx4g -jar snpEff.jar  -v  maizeZmB73_AGPV4  /XXX/vcf_file/Maize100_Chr"$i".vcf_maf5_miss20  -no-downstream -no-intergenic -no-upstream > /XXX/SNP_anno/Maize100_Chr"$i".vcf_maf5_miss20_anno_genic &
```

**For soybean:**

```
$	for i in {1..20}; do java -Xmx4g -jar snpEff.jar  -v  gmaxW82_v2  /XXX/vcf_file/Gm302_Chr"$i".vcf_maf5_miss20  -no-downstream -no-intergenic -no-upstream > /XXX/SNP_anno/Gm302_Chr"$i".vcf_maf5_miss20_anno_genic; done &
```

**(After this part, the procedures for maize and soybean are the same, and only need to change the input file, output file, and some database)**


**2. Divide genome-wide snps into genic and nongenic regions**

```
$	perl Extract_Maize100_genic_nongenic_SNPs.pl
```

**3. Divide nongenic SNPs into the gene-proximal region(promoter 1Kb and promoter 5Kb regions**

```
$ perl Divide_Maize100_nongenic_region_SNPs.pl
```	
	
**4. Calculate genome-wide base composition from genic, nongenic SNPs and the some other categories separately**

```
$ perl Cal_basecontent_from_chromosome_hmp_new_for_genic_nongenic_analysis_maf5miss20.pl
```

**5. Calculate base composition across chromosome segement for different different cateogories (Genic, nongenic., etc) and at different subtypes**

```
$ perl Cal_basecontent_across_chromosome_segment_for_genic_nongenic_analysis_all_category_maf5miss20.pl 1 10 &
```

**6. Random sample the same number of snps from both genic and nongenic to check the snp_effect**

```
$ perl Cal_basecontent_for_genic_nongenic_analysis_check_snp_number_effect.pl
``` 
	
**7. Calculate the MAF distribution**

```
$ Get_maf_for_genic_non_genic.pl
```

**8. Calculate base composition at nongenic and genic conditional on pericentromeric, and chromosome separately**

```
$  perl Cal_basecontent_for_genic_nongenic_conditional_on_peri_nonperi_analysis_maf5miss20.pl
```

**9. Calculate base composition at selective sweep and non-selective sweep regions**

```
$ perl Cal_basecontent_for_sweep_unsweep_analysis_maf5miss20.pl &
```