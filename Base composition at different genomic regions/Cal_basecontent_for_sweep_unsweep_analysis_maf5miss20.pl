#!/usr/local/bin/perl -w
###This code filters the chromosomm_hmp_file to certain maf and missing rate. 
###And it also calculate the ACGT base content for each strain across 5 chromosome. 

use strict;

#####################Specifying the parameters
my $group_info = "Maize100";
my ($ch_s, $ch_e) = (1, 10);
my $sweep_snpn = 635240;
my $unsweep_snpn =  8217438;
my $missing_percent = 20;
my $maf_percent = 5;


#####################Specifying the direcotry
my $unsweep_sample_percentage = sprintf "%.6f", $sweep_snpn/$unsweep_snpn;

my $dir = '/XXX/';

my $input_snp_file_dir = $dir.$group_info.'/geno/';
my $sweep_file_dir = $dir.$group_info.'/Selective_sweep_analysis/';

my $id_line_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_1st_id_line';
my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 

my $strain_group_infor_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_group_infor';
my $strain_group_hashref = Get_group_infor($strain_group_infor_file);

my @category_array	= ("All", "sweep", "unsweep", "sampled_unsweep" );

my @bases = qw /A C G T/;
my (%strain_base_count, %strain_snp_count);

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	
	my $input_chrsnp_file = $input_snp_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 	
	next unless -e $input_chrsnp_file;
	my $sweep_snp_file = $sweep_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent.'_sweep_infor';
	next unless -e $sweep_snp_file;
	
	my $snp_info_hashref = Get_snp_infor($sweep_snp_file);
	
	open (INPUT, $input_chrsnp_file ) || die;
#	my $k = 0;
		while (<INPUT>) {
			chomp;
			next if /alleles/;
			my @t = split /\t/;
#			$k++;
#			next if $k > 10000;
			my $bp = $t[3];
			my $flag1 = 0;
			my $category = $$snp_info_hashref{$bp};
##This one is to sample certain number of unsweep snps			
			if ($$snp_info_hashref{$bp} eq "unsweep" && rand() <= $unsweep_sample_percentage){
#				next unless (rand() <= $unsweep_sample_percentage);
				$flag1 = 1;
			}

			for (my $i = 11; $i <= $#t; $i ++) {
				next unless $t[$i] =~ /[ACGT]/;
				my $strain = $$strain_id_arrayref[$i];
				$strain_base_count{$category_array[0]}{$strain}{$t[$i]} ++;
				$strain_snp_count{$category_array[0]}{$strain} ++;
				$strain_base_count{$category}{$strain}{$t[$i]} ++;
				$strain_snp_count{$category}{$strain} ++;
				if ($flag1 == 1) {
					$strain_base_count{$category_array[3]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[3]}{$strain} ++;
					}
				}
			}
		close INPUT;	
#		close Ch;
	}
		
my $output_dir = $sweep_file_dir.'BCS/';
mkdir $output_dir unless -e $output_dir; 			
my $output_file = $output_dir.$group_info.'_sweep_nonsweep_ATCG_maf'.$maf_percent.'_miss'.$missing_percent;
open (OUT, '>'.$output_file) || die;
print OUT "Category\tStrain\tGroup\tA\tC\tG\tT\tTotal\n";

foreach my $categ (@category_array) {
	
	for (my $i = 11; $i < @$strain_id_arrayref; $i ++)  {
		my $strain = $$strain_id_arrayref[$i];
		my $group = $$strain_group_hashref{$strain};
		print OUT $categ."\t".$strain."\t".$group;
		my $total_num = exists $strain_snp_count{$categ}{$strain} ? $strain_snp_count{$categ}{$strain} : 1;
		foreach my $b (@bases) {
			my $base_count = exists $strain_base_count{$categ}{$strain}{$b} ? $strain_base_count{$categ}{$strain}{$b} : 0;
			my $percent = sprintf "%.5f", $base_count / $total_num;
			print OUT "\t".$percent;
		}
		print OUT "\t".$total_num."\n";
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
	my ($sweep_f) = @_;
	my %snp_infor;
	open (F1, $sweep_f) || die;
		while (<F1>) {
		chomp;
		next if /Allele/;
		my @t = split /\t/;
		$snp_infor{$t[3]} = $t[4];  
		}
		close F1;
		return \%snp_infor;
	}