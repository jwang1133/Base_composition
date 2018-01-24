#!/usr/local/bin/perl -w

use strict;

##Specify the input file directory
my $ori_dir = '/XXX/';
##Speficy the output file directory
my $hmp_dir = '/XXX/';

my $group = 'Gm';
my $group_number = 302;
my $missing_percent = 100;
my $maf_percent = 5;

my $missing_rate = $missing_percent/100;
my $maf = $maf_percent/100;
my $group_info = $group.$group_number;

##User input for the specific chromosomes need to loop through
my ($ch_s, $ch_e) = @ARGV;
for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $hmp_file = $hmp_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 
	my $chch = $ch < 10 ? '0'.$ch : $ch;
	my $ori_ch_file = $ori_dir.'Gm'.$chch.'.snp.genotype';
	Format_to_hmp($ori_ch_file, $hmp_file, $chch, $group_number, $missing_rate, $maf);
	}

sub Format_to_hmp {
	my ($f, $o, $chro, $group_n, $missing_r, $maf_r) = @_;
	open (F, $f) || die;
	open (Hmp, '>'.$o) || die;
	my @hmp_pre = qw/rs alleles chrom pos strand assembly center protLSID assayLSID panel/;
	print Hmp '';
	print Hmp $_."\t" foreach (@hmp_pre);
	print Hmp 'QC_code';
	print Hmp "\tGm_".$_ foreach (1..302);
	print Hmp "\n";

	while (<F>) {
		chomp;
		my @t = split /\s+/;
		if ($t[0] eq 'Chromosome') {
			}
			
			else {
				my $ch = $t[0];
				my $bp = $t[1];	
				my @bp_allele;
			
				for (my $i = 2; $i <= $#t; $i ++)	{
				my $g;
				if ($t[$i] =~ /[ACGT]/) {$g = $t[$i];}
				else {$g = 'N';}
				push @bp_allele, $g;	
					}	
				my $total=0;  
				my %hash;
				foreach my $a (@bp_allele) {
					next if $a =~ /N/;
					$total ++;
					$hash{$a} ++;
				}
					
				next unless $total >= (1-$missing_r) * $group_n;
				my @sorted_allele	= sort { $hash{$a} <=> $hash{$b} } keys %hash;
				##control if there is three homozygous allele at one loci, ignore it.
				next unless $#sorted_allele == 1;  
				my $maf = sprintf "%.3f", $hash{$sorted_allele[0]} / $total;
				next unless $maf >= $maf_r;
				my $ref_allele = $sorted_allele[1]; 
				##Here $alt_allele is minor allele.
				my $alt_allele = $sorted_allele[0]; 
				my $ref_alt_allele = join('/', $ref_allele,$alt_allele);			
				my $pos_x = sprintf "%08d", $bp;
				my $rs_id = 'rs'.$chro.$pos_x;
				print Hmp $rs_id."\t".$ref_alt_allele."\t".$ch."\t".$bp."\t+\tGm302\tNA\t$maf\t$total\tNA\tNA";
				print Hmp "\t".$_ foreach (@bp_allele);
				print Hmp "\n";		
			}	
		}
		close F;
	}
