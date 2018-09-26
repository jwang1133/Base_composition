#!/usr/local/bin/perl -w
###This code filters the chromosomm_hmp_file to certain maf and missing rate. 
###And it also calculate the ACGT base content for each strain across 5 chromosome. 

use strict;
## merge segment hmps into chr
##those are the parameter need to change

#######################################
my $group = 'Maize';
my $group_number=100;
my $missing_percent = 20;
my $maf_percent = 5;
my $group_info = $group.$group_number;

#####################################
#my @subs_type = qw/AC AG AT CG CT GT/;

my $sub = 'AC';
#my $sub = 'AG';
#my $sub = 'AT';
#my $sub = 'CG';
#my $sub = 'CT';
#my $sub = 'GT';

#######################################
my $input_file_dir = '/XXX/';
my $output_dir_pheno = '/XXX/';

my $id_line_file = $input_file_dir.$group_info.'_1st_id_line';
my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 

my $total_snp_file = 'Maize100_ATCG_maf5_miss20';
my $snp_numer_hashref = Get_snp_number($total_snp_file); 

my @bases = qw /A C G T/;
my %strain_base_count;


for (my $ch = 1; $ch <= 10; $ch ++) {

	my $input_chrsnp_file = $input_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 
	next unless -e $input_chrsnp_file;	
	open (INPUT, $input_chrsnp_file ) || die;
	
		while (<INPUT>) {
			chomp;
			my $line = $_;
			my @t = split /\t/, $line;
			next if $t[1] eq 'alleles';
			my ($a1, $a2) = split /\//, $t[1];
		  my $a1a2 = join '', sort($a1, $a2);	
		  next unless $a1a2 eq $sub;	
			for (my $i = 11; $i <= $#t; $i ++) {
				next unless $t[$i] =~ /[ACGT]/;
				my $strain = $$strain_id_arrayref[$i];
				$strain_base_count{$strain}{$t[$i]} ++;
				}
			}
		close INPUT;	
	}
		
my $output_pheno = $output_dir_pheno.$group_info.'_Subgroup_'.$sub.'_ATCG_maf'.$maf_percent.'_miss'.$missing_percent;

open (OUT, '>'.$output_pheno) || die;
print OUT "Strain\tA\tC\tG\tT\n";
for (my $i = 11; $i < @$strain_id_arrayref; $i ++)  {
	my $strain = $$strain_id_arrayref[$i];
	print OUT $strain;
	my $total_num = exists $$snp_numer_hashref{$strain} ? $$snp_numer_hashref{$strain} : 1;
	foreach my $b (@bases) {
		my $base_count = exists $strain_base_count{$strain}{$b} ? $strain_base_count{$strain}{$b} : 0;
		my $percent = sprintf "%.5f", $base_count / $total_num;
		print OUT "\t".$percent;
		}
	print OUT "\n";
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
	
	
sub Get_snp_number {
	my $file = shift;
	open (FILE, $file ) || die;
	my %hash;
	while (<FILE>) {
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'Strain';
		$hash{$t[0]} = $t[5];
		}
		close FILE;
		return \%hash;
	}	
	
	
