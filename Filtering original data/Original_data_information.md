## This file describes the original data source and how the genotype data is converted to the coordinates of the most advanced version of reference genome##

**1. Maize original SNP data**

-Maize original dataset can be downloaded from website [http://www.panzea.org/genotypes](http://www.panzea.org/genotypes) 

After open the above link, click on the link for Maize HapMapV2 genotypes and then download the SNP file in hapmap format 

**2.Soybean original SNP data**

-Soybean original dataset can be downloaded from website [http://figshare.com/articles/Soybean_resequencing_project/1176133](http://figshare.com/articles/Soybean_resequencing_project/1176133)

After open the link download the file named SNP-1.zip, SNP-2.zip and SNP-3.zip

**3. Update the maize genotype data in B73 RefV2 to B73 RefV4**



-  Download the chain file from website [http://ftp.gramene.org/release-57/assembly_chain/zea_mays/AGPv2_to_AGPv4.chain.gz](http://ftp.gramene.org/release-57/assembly_chain/zea_mays/AGPv2_to_AGPv4.chain.gz)


-  Prepare the bed file for the genotype with the script named **prepare\_maize\_bed\_file\_new.pl**
 

- Convert the genotype coordinate in B73 RefV2 to B73 RefV4 with **CrossMap**
```
 $ for i in {1..10}; do python /XXX/CrossMap.py  bed  /XXX/AGPv2_to_AGPv4.chain.gz  /XXX/maizeHapMapV2_B73RefGenV2_201203028_chr"$i"_v2_coordinate  >  /XXX/maizeHapMapV2_B73RefGenV2_201203028_chr"$i"_v2_v4_coordinate ; done &
```
- Run the the following script to convert the hmp file in v2 to hmp file in v4
``` 
$ perl  convert_maize_v2_to_v4_original.pl
```
- Concatenate all the converted chromosome file to be one file, and then reprint each of the converted hmp file to make sure that SNPs on the same chromsome are together in one file
```
$ for i in {1..10}; do cat maizeHapMapV2_B73RefGenV2_201203028_chr"$i"_v4.hmp.txt >> maizeHapMapV2_B73RefGenV2_201203028_genome_v4.hmp.txt; done &
```
```
$ awk '$3==1' maizeHapMapV2_B73RefGenV2_201203028_genome_v4.hmp.txt | sort -k4,4n > maizeHapMapV2_B73RefGenV2_201203028_chr1_v4_sorted.hmp.txt &
```
**4. Update the soybean genotype data glycine\_max v1.0 to v2.0** 

- Download the chain file from [http://ftp.gramene.org/release-57/assembly_chain/glycine_max/V1.0_to_Glycine_max_v2.0.chain.gz](http://ftp.gramene.org/release-57/assembly_chain/glycine_max/V1.0_to_Glycine_max_v2.0.chain.gz)

- The step for preparing the bed file and convert the coordinates with CrossMap function are similar as that for maize
- Then Run the following script to convert the soybean hmp file in v1 to hmp file in v2
```
$ perl  convert_soybean_v1_to_v2_original.pl
```
- The **concatenate and reprint** procedures are the same as that in maize