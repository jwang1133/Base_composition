#!/usr/local/bin/perl -w
##This script aims to convert the hmp file in W82 RefV1 to RefV2 version. 

use strict;
my $dir = '/XXX/';
my $v1_dir = $dir.'original/'; 
my $v2_dir = $dir.'original/original_in_v2/';
mkdir $v2_dir unless -e $v2_dir;
my $v1_v2_dir = $dir.'v1_v2/';

my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

##This step read in all the v2_v4 coordinate file to a hash%, so that we can have all the v4 chromosome information
##And print out the first line for each output file in v4
my %hash;
for (my $ch = $ch_s; $ch <= $ch_e; $ch ++){
	my $v1_v2_file = $v1_v2_dir.'Gm'.$ch.'.snp.genotype_v1_v2_coordinate';
	open (F, $v1_v2_file) || die;
	while (<F>) {
		chomp;
		my $line = $_;
		next if $line =~ /Fail/;
		my @t = split /\t/;
		next unless $t[4] =~ /^\d+?$/; ##check whether the SNP mactchs to any of the 10 chromsomes
		$hash{$t[0]."\t".$t[1]} =$t[4]."\t".$t[5];
		}
	close F;	
	my $chch = $ch < 10 ? '0'.$ch : $ch;
	my $ch_v1_file = $v1_dir.'Gm'.$chch.'.snp.genotype';
	my $ch_v2_file = $v2_dir.'Gm'.$ch.'.snp.genotype_v2';
	open (V2, '>'.$ch_v2_file) || die;
	open (V1, $ch_v1_file) || die;
	while (<V1>) {
		chomp;
		my $line = $_;
		next if $line =~ /alleles/;
		my @t = split /\s+/;
		next unless exists $hash{$t[0]."\t".$t[1]};
		my $new_snp_pos = $hash{$t[0]."\t".$t[1]};
#		my @new_snp_pos = split(/\t/,$new_snp_pos);
#		my $chrom = $new_snp_pos[0];	
		print V2  $new_snp_pos."\t";
		print V2  $_."\t" foreach (@t[2..303]);
		print V2 "\n";
				}
	close V1;
	close V2;
	}
	

#sub Parse_crossmap_file {
#	my ($f) = @_;
#	my (%hash);
#	open (F, $f) || die;
#	while (<F>) {
#		chomp;
#		my $line = $_;
#		next if $line =~ /Fail/;
#		my @t = split /\t/;
#		next unless $t[0] eq $t[4];
#		$hash{$t[1]} = $t[5];
#		}
#	close F;
#	return \%hash;	
#	}	