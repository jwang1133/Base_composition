#!/usr/local/bin/perl -w

use strict;

my $group = 'Gm';
my $group_number=302;
my $missing_percent = 20;
my $maf_percent = 5;

my $missing_rate = $missing_percent/100;
my $maf = $maf_percent/100;
my $group_info = $group.$group_number;


##Specify the input and output file directory
my $input_file_dir = '/XXX/';
my $output_dir_pheno = '/XXX/';


my $id_line_file = $input_file_dir.'Gm302_1st_id_line';

my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 
my @bases = qw /A C G T/;
my (%strain_base_count, %strain_snp_count);

for (my $ch = 1; $ch <= 20; $ch ++) {

	my $input_chrsnp_file = $input_file_dir.'Gm302_Chr'.$ch.'.hmp_maf5_miss20';
	next unless -e $input_chrsnp_file;	
	open (INPUT, $input_chrsnp_file ) || die;
	my %bp_snp_info;
		while (<INPUT>) {
			chomp;
			my $line = $_;
			next if $line =~ /alleles/; # *If it is the first row in the file, Jump to next line
			my @t = split /\t/, $line;		
			for (my $i = 11; $i <= $#t; $i ++) {
				next unless $t[$i] =~ /[ACGT]/;
				my $strain = $$strain_id_arrayref[$i];
				$strain_base_count{$strain}{$t[$i]} ++;
				$strain_snp_count{$strain} ++;
				}
			}
		close INPUT;	
	}
		
my $output_pheno = $output_dir_pheno.$group_info.'_ATCG_maf'.$maf_percent.'_miss'.$missing_percent;

open (OUT, '>'.$output_pheno) || die;
print OUT "Strain\tA\tC\tG\tT\tTotal\n";

for (my $i = 11; $i < @$strain_id_arrayref; $i ++)  {
	my $strain = $$strain_id_arrayref[$i];
	print OUT $strain;
	my $total_num = exists $strain_snp_count{$strain} ? $strain_snp_count{$strain} : 1;
	foreach my $b (@bases) {
		my $base_count = exists $strain_base_count{$strain}{$b} ? $strain_base_count{$strain}{$b} : 0;
		my $percent = sprintf "%.5f", $base_count / $total_num;
		print OUT "\t".$percent;
		}
	print OUT "\t".$total_num."\n";
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
	
	
