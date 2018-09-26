#!/usr/local/bin/perl -w
use strict;

my $group_info = 'Maize100';
my $missing_percent = 20;
my $maf_percent = 5;
my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

my $dir = '/hdd/jinyuw/BCS/';


my $input_snp_file_dir = $dir.$group_info.'/geno/';
my $genic_non_genic_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';

my $output_dir = $genic_non_genic_file_dir.'MAF/';
mkdir $output_dir unless -e $output_dir; 			
my $output_file = $output_dir.$group_info.'_MAF_genic_nongenic_maf'.$maf_percent.'_miss'.$missing_percent.'_Chr'.$ch_s.'_'.$ch_e;

open (OUT, '>'.$output_file) || die;
print OUT "Category\tChromosome\tPosition\tAllele\tMissing_rate\tMAF\n";


	for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	
		my $input_chrsnp_file = $input_snp_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 	
		next unless -e $input_chrsnp_file;
		my $nongenic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic';
		next unless -e $nongenic_snp_file;
		my $snp_info_hashref = Get_snp_infor($nongenic_snp_file);
		
		open (INPUT, $input_chrsnp_file ) || die;
#		my $k = 0;
			while (<INPUT>) {
				chomp;
				next if /rs#/;
		  	my @t = split /\t/;
#				$k++;
#				last if $k > 10000;
				my ($allele,$bp, $missingr, $maf) = ($t[1], $t[3], $t[9], $t[10]);
				if (exists $$snp_info_hashref{$bp}){
						print OUT "Nongeic\t".$ch."\t".$bp."\t".$allele."\t".$missingr."\t".$maf."\n";
					}
				else {
						print OUT "Genic\t".$ch."\t".$bp."\t".$allele."\t".$missingr."\t".$maf."\n";
					} 
			}
	close INPUT;	
}	





	sub Get_snp_infor {
	my ($nongenic_f) = @_;
	my %snp_infor;
		open (F2, $nongenic_f ) || die;
		while (<F2>) {
		chomp;
		next if /Allele/;
		my @t = split /\t/;
		$snp_infor{$t[3]} = $t[4];  
		}
		close F2;
		return \%snp_infor;
	}