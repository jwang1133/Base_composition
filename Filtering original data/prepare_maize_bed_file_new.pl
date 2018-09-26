#!/usr/local/bin/perl -w
use strict;

my $input_file_dir = '/XXX/';
my $output_file_dir = '/XXX/';
mkdir $output_file_dir unless -e $output_file_dir;
my ($ch_s, $ch_e) = (1, 10);

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $input_file = $input_file_dir.'maizeHapMapV2_B73RefGenV2_201203028_chr'.$ch.'.hmp.txt';
	next unless -e $input_file;
	
	my $output_file = $output_file_dir.'maizeHapMapV2_B73RefGenV2_201203028_chr'.$ch.'_v2_coordinate';

	open (F, $input_file)||die;
	open (O, '>'.$output_file) || die;
	print O "Chromosome"."\t"."Pos_s"."\t"."Pos_e"."\n";
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'rs#';
		print O $t[2]."\t".$t[3]."\t".$t[3]."\n";
		}
	close F;
}