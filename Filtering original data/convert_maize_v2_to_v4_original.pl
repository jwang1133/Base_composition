#!/usr/local/bin/perl -w
##This script aims to convert the hmp file in B73 RefV2 to RefV4 version. 

use strict;
my $dir = '/XXX/';
my $v2_dir = $dir.'original/'; 
my $v4_dir = $dir.'original/original_in_v4/';
mkdir $v4_dir unless -e $v4_dir;
my $v2_v4_dir = $dir.'v2_v4/';

my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

##This step read in all the v2_v4 coordinate file to a hash%, so that we can have all the v4 chromosome information
##And print out the first line for each output file in v4
my %hash;
for (my $ch = $ch_s; $ch <= $ch_e; $ch ++){
	my $v2_v4_file = $v2_v4_dir.'maizeHapMapV2_B73RefGenV2_201203028_chr'.$ch.'_v2_v4_coordinate';
	open (F, $v2_v4_file) || die;
	while (<F>) {
		chomp;
		my $line = $_;
		next if $line =~ /Fail/;
		my @t = split /\t/;
		next unless $t[4] =~ /^\d+?$/; ##check whether the SNP mactchs to any of the 10 chromsomes
		$hash{$t[0]."\t".$t[1]} =$t[4]."\t".$t[5];
		}
	close F;	
	
	my $ch_v2_file = $v2_dir.'maizeHapMapV2_B73RefGenV2_201203028_chr'.$ch.'.hmp.txt';
	my $ch_v4_file = $v4_dir.'maizeHapMapV2_B73RefGenV2_201203028_chr'.$ch.'_v4.hmp.txt';
	open (V4, '>'.$ch_v4_file) || die;
	open (V2, $ch_v2_file) || die;
	while (<V2>) {
		chomp;
		my $line = $_;
		next if $line =~ /alleles/;
		my @t = split /\t/;
		$t[5] = 'RefGenV4';
		next unless exists $hash{$t[2]."\t".$t[3]};
		my $new_snp_pos = $hash{$t[2]."\t".$t[3]};
#		my @new_snp_pos = split(/\t/,$new_snp_pos);
#		my $chrom = $new_snp_pos[0];
		print V4  $_."\t" foreach (@t[0..1]);
		print V4  $new_snp_pos;
		print V4 "\t".$_ foreach (@t[4..114]);
		print V4 "\n";
				}
	close V2;
	close V4;
	}

