#!/usr/local/bin/perl -w
###This code filters the chromosomm_hmp_file to certain maf and missing rate. 
###And it also calculate the ACGT base content for each strain across 5 chromosome. 
use strict;


###############Specifying the parameters
my $group_info = 'Maize100';
my $missing_percent = 20;
my $maf_percent = 5;

##Plan to sample about 90% of the snps.
my $sampled_snp_number = 1558518;
my $genic_snp_number = 1731687;
my $nongeic_snp_number = 7120981;
#my $genic_r = 0.1956;
#my $nongenic_r = 0.8044;

###############Specifying the segment and step infomration
#my $sample_percentage = sprintf "%.6f", $genic_r / $nongenic_r;
my $sample_genic_percentage = sprintf "%.6f", $sampled_snp_number / $genic_snp_number;
my $sample_nongenic_percentage = sprintf "%.6f", $sampled_snp_number / $nongeic_snp_number;

my ($round_s, $round_e ) = (1, 100);

###############Directory information 
my $dir = '/XXX/';

my $input_snp_file_dir = $dir.$group_info.'/geno/';
my $genic_non_genic_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';


my $id_line_file = $input_snp_file_dir.$group_info.'_1st_id_line';
my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 

my $strain_group_infor_file = $input_snp_file_dir.$group_info.'_group_infor';
my $strain_group_hashref = Get_group_infor($strain_group_infor_file);

my @category_array = ("Genic", "Nongenic");

my $output_dir = $genic_non_genic_file_dir.'BCS/';
mkdir $output_dir unless -e $output_dir; 			
my $output_file = $output_dir.$group_info.'_genic_nongenic_ATCG_check_snp_number_effect_maf'.$maf_percent.'_miss'.$missing_percent.'_iteration'.$round_e;
open (OUT, '>'.$output_file) || die;
print OUT "Iteration\tCategory\tStrain\tGroup\tA\tC\tG\tT\tTotal\n";


my @bases = qw /A C G T/;
my (%strain_base_count, %strain_snp_count);

for (my $round = $round_s; $round <= $round_e; $round ++ ){
	for (my $ch = 1; $ch <= 10; $ch ++) {
	
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
#				next if $k > 100;
				my $bp = $t[3];
				my ($flag1, $flag2) = qw/0 0/;
#				next unless exists $$snp_info_hashref{$bp};
#				print $$snp_info_hashref{$bp}."\n";
				if ($$snp_info_hashref{$bp}){
					next unless (rand()<= $sample_nongenic_percentage);
					$flag1 = 1;
				}
				else {
					next unless (rand()<= $sample_genic_percentage);
					$flag2 = 1;
				}
				for (my $i = 11; $i <= $#t; $i ++) {
					next unless $t[$i] =~ /[ACGT]/;
					my $strain = $$strain_id_arrayref[$i];
					if ($flag1 == 1){
						$strain_base_count{$round}{$category_array[1]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$round}{$category_array[1]}{$strain} ++;
						}
					if ($flag2 == 1){
						$strain_base_count{$round}{$category_array[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$round}{$category_array[0]}{$strain} ++;
						}
				}
			}
			close INPUT;	
	}
	
	foreach my $categ (@category_array) {
		for (my $i = 11; $i < @$strain_id_arrayref; $i ++)  {
			my $strain = $$strain_id_arrayref[$i];
			my $group = $$strain_group_hashref{$strain};
			print OUT $round."\t".$categ."\t".$strain."\t".$group;
			my $total_num = exists $strain_snp_count{$round}{$categ}{$strain} ? $strain_snp_count{$round}{$categ}{$strain} : 1;
			foreach my $b (@bases) {
				my $base_count = exists $strain_base_count{$round}{$categ}{$strain}{$b} ? $strain_base_count{$round}{$categ}{$strain}{$b} : 0;
				my $percent = sprintf "%.5f", $base_count / $total_num;
				print OUT "\t".$percent;
			}
			print OUT "\t".$total_num."\n";
		}	
	}
}		











sub Get_strain_id {
	my $file = shift;
	open (FILE, $file ) || die;
	my (@array, $line);
	while (<FILE>) {
		chomp;
			$line = $_;
		  my @t = split /\t/, $line;
		  @array = @t;  
		}
		close FILE;
		return (\@array, $line);
	}	
	
sub Get_group_infor {
	my $file = shift;
	open (FILE, $file ) || die;
	my %hash;
	while (<FILE>) {
		chomp;
		next if /Group/;
		my @t = split /\t/;
		$hash{$t[0]} = $t[1];  
		}
		close FILE;
		return (\%hash);
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