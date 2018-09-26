#!/usr/local/bin/perl -w

use strict;

###############Specifying the parameters
my $group_info = "Maize100";
my ($ch_s, $ch_e) = (1, 10);
my $missing_percent = 20;
my $maf_percent = 5;

###############Specifying the directory

my $dir = '/XXX/';

my $cent_infor_hashref = Parse_cent_infor($dir.$group_info.'/Genome_db/'.$group_info.'_pericent_cent_infor');

my $input_snp_file_dir = $dir.$group_info.'/geno/';
my $genic_non_genic_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';

my $id_line_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_1st_id_line';
my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 

my $strain_group_infor_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_group_infor';
my $strain_group_hashref = Get_group_infor($strain_group_infor_file);



my @category_array	= ("All", "peri", "nonperi", "peri_genic", "peri_nongenic", "nonperi_genic", "nonperi_nongenic" );

my @bases = qw /A C G T/;
my (%strain_base_count, %strain_snp_count);

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
##Get pericetnromere information	
	my $cent_infor = $$cent_infor_hashref{$ch};
	my ($pericent_s, $pericent_e) = split /\t/, $cent_infor;
	
	my $input_chrsnp_file = $input_snp_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 	
	next unless -e $input_chrsnp_file;
	my $nongenic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic';
	next unless -e $nongenic_snp_file;
	my $snp_info_hashref = Get_snp_infor($nongenic_snp_file);
	
	
	open (INPUT, $input_chrsnp_file ) || die;
#	my $k = 0;
		while (<INPUT>) {
			chomp;
			next if /alleles/;
			my @t = split /\t/;
#			$k++;
#			next if $k > 10000;
			my $bp = $t[3];
			my ($flag1, $flag2, $flag3, $flag4, $flag5, $flag6 ) = qw/0 0 0 0 0 0/;
			my $category = $$snp_info_hashref{$bp};
##This part is to check which category the snp belongs to 
			if ($bp >= $pericent_s && $bp <= $pericent_e ){
				$flag1 = 1;
			}
			else {
				$flag2 = 1;
			}
			
			if ($$snp_info_hashref{$bp}) {
				if ($flag1 == 1){
					$flag4 = 1;
				}
				if ($flag2 == 1){
					$flag6 = 1;
				}
			}
			else {
				if ($flag1 == 1){
					$flag3 = 1;
				}
				if ($flag2 == 1){
					$flag5 = 1;
				}
			}
			
			for (my $i = 11; $i <= $#t; $i ++) {
				next unless $t[$i] =~ /[ACGT]/;
				my $strain = $$strain_id_arrayref[$i];
				$strain_base_count{$category_array[0]}{$strain}{$t[$i]} ++;
				$strain_snp_count{$category_array[0]}{$strain} ++;			
				if ($flag1 == 1) {
					$strain_base_count{$category_array[1]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[1]}{$strain} ++;
				}
				if ($flag2 == 1) {
					$strain_base_count{$category_array[2]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[2]}{$strain} ++;
				}
				if ($flag3 == 1) {
					$strain_base_count{$category_array[3]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[3]}{$strain} ++;
				}
				if ($flag4 == 1) {
					$strain_base_count{$category_array[4]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[4]}{$strain} ++;
				}
				if ($flag5 == 1) {
					$strain_base_count{$category_array[5]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[5]}{$strain} ++;
				}
				if ($flag6 == 1) {
					$strain_base_count{$category_array[6]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[6]}{$strain} ++;
				}
			}
		}
		close INPUT;	
#		close Ch;
	}
		
my $output_dir = $genic_non_genic_file_dir.'BCS/';
mkdir $output_dir unless -e $output_dir; 			
my $output_file = $output_dir.$group_info.'_peri_nonperi_genic_nongenic_ATCG_maf'.$maf_percent.'_miss'.$missing_percent;
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
	
	
sub Parse_cent_infor {
	my ($file) = @_;
	my %hash;
	open (FILE, $file) || die;
	while (<FILE>) {
		chomp;
		next if /Chromosome/;
		my @t = split /\t/;
		$hash{$t[0]} = $t[1]."\t".$t[2];
		}
		close FILE;
		return \%hash;
	}	