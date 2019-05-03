#!/usr/local/bin/perl -w
###This code filters the chromosomm_hmp_file to certain maf and missing rate. 
###And it also calculate the ACGT base content for each strain across 5 chromosome. 

use strict;
## merge segment hmps into chr
##those are the parameter need to change

my $group = 'Maize';
my $group_number = 100;
##For soybean
#my $group = 'Gm';
#my $group_number = 302;
my $missing_percent = 20;
my $maf_percent = 5;
###################################
my $group_info = $group.$group_number;

my $dir = '/XXX/';


my $input_snp_file_dir = $dir.$group_info.'/geno/';
my $genic_non_genic_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';

my $id_line_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_1st_id_line';
my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 

my $strain_group_infor_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_group_infor';
my $strain_group_hashref = Get_group_infor($strain_group_infor_file);

my @category_array	= ("All", "Genic", "Nongenic","Nonsyn","UTRs","Og","Promoter5k","In", "Syn", "Nyn", "3U", "5U", "Sg", "Rg","Non",  "P5", "P1");

my @Genic_array = ("In", "Syn", "Nyn", "3U", "5U", "Sg", "Rg" );
my @Nongenic_array = ("Non", "P5", "P1"); 
####Different ways to classify the Nonsynonymous SNPs: "Nonsyn"#####
my @Nsy_array = ("Nyn", "Sg", "Rg" );
#############################"UTRs"#####################
my @UTRs_array = ( "3U", "5U");
#############################"Og"#####################
my @Othergenic_array = ("Sg", "Rg" );

my @Promoter5K_array = ("P5", "P1");


my @bases = qw /A C G T/;
my (%strain_base_count, %strain_snp_count);

for (my $ch = 1; $ch <= 10; $ch ++) {
	
	my $input_chrsnp_file = $input_snp_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 	
	next unless -e $input_chrsnp_file;
	my $genic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_genic';
	next unless -e $genic_snp_file;
	my $nongenic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic_splitted';
	next unless -e $nongenic_snp_file;
	
	my $snp_info_hashref = Get_snp_infor($genic_snp_file, $nongenic_snp_file);
	
	open (INPUT, $input_chrsnp_file ) || die;
#	my $k = 0;
		while (<INPUT>) {
			chomp;
			next if /alleles/;
			my @t = split /\t/;
#			$k++;
#			next if $k > 10000;
			my $bp = $t[3];
			my ($flag1, $flag2, $flag3, $flag4, $flag5, $flag6) = qw/0 0 0 0 0 0/;
			my $category = $$snp_info_hashref{$bp};
#			print $category."\n";
			if (grep {$_ eq $category} @Genic_array) {
 				$flag1 = 1;
			}
			if (grep {$_ eq $category} @Nongenic_array) {
 				$flag2 = 1;
			}
			if (grep {$_ eq $category} @Nsy_array) {
 				$flag3 = 1;
			}
			
			if (grep {$_ eq $category} @UTRs_array) {
 				$flag4 = 1;
			}		
			if (grep {$_ eq $category} @Othergenic_array) {
 				$flag5 = 1;
			}	
			if (grep {$_ eq $category} @Promoter5K_array) {
 				$flag6 = 1;
			}						
			for (my $i = 11; $i <= $#t; $i ++) {
				next unless $t[$i] =~ /[ACGT]/;
				my $strain = $$strain_id_arrayref[$i];
				$strain_base_count{$category_array[0]}{$strain}{$t[$i]} ++;
				$strain_snp_count{$category_array[0]}{$strain} ++;
				$strain_base_count{$category}{$strain}{$t[$i]} ++;
				$strain_snp_count{$category}{$strain} ++;
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
my $output_file = $output_dir.$group_info.'_genic_nongenic_ATCG_maf'.$maf_percent.'_miss'.$missing_percent;
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
	my ($genic_f, $nongenic_f) = @_;
	my %snp_infor;
		open (F1, $genic_f ) || die;
		while (<F1>) {
		chomp;
		next if /Allele/;
		my @t = split /\t/;
		$snp_infor{$t[3]} = $t[4];  
		}
		close F1;
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