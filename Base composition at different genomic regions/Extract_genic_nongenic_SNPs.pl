#!/usr/local/bin/perl -w

##This file extract the genic and non_genic snps

use strict;

#############Specifying the parameters########
my $group = 'Maize';
my $group_number=100;
my $missing_percent = 20;
my $maf_percent = 5;
###################################
my $group_info = $group.$group_number;
my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

#############Specifying the directory##########

my $dir = '/XXX/';
my $input_file_dir = $dir.$group_info.'/SNP_anno/';
my $output_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';
mkdir $output_file_dir unless -e $output_file_dir;



for (my $ch = $ch_s; $ch <= $ch_e; $ch ++){
#	my $chch = sprintf "%02d", $ch;
	my $input_file = $input_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_anno_genic';
	next unless -e $input_file;
	my $output_file1 = $output_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_genic';
	my $output_file2 = $output_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic';
	
	open (O1, '>'.$output_file1) || die;
	print O1 "SNP\tAllele\tChromosome\tPos\tCategory\n";
	
	open (O2, '>'.$output_file2) || die;
	print O2 "SNP\tAllele\tChromosome\tPos\tCategory\n";
	
	open (In, $input_file) || die;
	while (<In>) {
		chomp;
		next if /#/;
#		next if $k > 1000;
#		$k++;
		my @t = split /\t/;
		my $allele_info = join('/', $t[3], $t[4]);
		my @snp_anno = split(';', $t[7]);
		if ($#snp_anno == 0){
			print O2 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."Non"."\n";
			}
		else {
			my @snp_infor = split(/\|/, $snp_anno[1]);
#			print $snp_infor[1]."\n";
			if ($snp_infor[1] =~ /intron_variant/){
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."In"."\n";
				} elsif ($snp_infor[1] =~ /synonymous_variant/){
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."Syn"."\n";
				}	elsif ($snp_infor[1] =~ /missense_variant/){
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."Nyn"."\n";			
			}		elsif ($snp_infor[1] =~ /3_prime_UTR/){
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."3U"."\n";			
			}		elsif ($snp_infor[1] =~ /5_prime_UTR/){
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."5U"."\n";			
			}		elsif ($snp_infor[1] =~ /stop_gained/){
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."Sg"."\n";			
			}		else {
					print O1 $t[2]."\t".$allele_info."\t".$ch."\t".$t[1]."\t"."Rg"."\n";			
			}
		}
	}
	close In;
}