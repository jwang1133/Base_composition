#!/usr/local/bin/perl -w

##This file extract the genic and non_genic snps

use strict;

#############Specifying the parameters########
my $group = 'Maize';
my $group_number = 100;
my $missing_percent = 20;
my $maf_percent = 5;
###################################
my $group_info = $group.$group_number;
my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

#############Specifying the directory########
my $dir = '/XXX/';
my $file_dir = $dir.$group_info.'/Genic_nongenic_analy/';
my $gene_file = $dir.$group_info.'/XXX/Zea_mays.AGPv4.39.chr.gff3';

my ($w1, $w2) = (1000, 5000);


for (my $ch = $ch_s; $ch <= $ch_e; $ch ++){
	my $input_file = $file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic';
	next unless -e $input_file;	
	my ($gene_pos_array1, $gene_pos_array2) = Get_gene_pos($gene_file, $ch);	
	my $output_file = $file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic_splitted';
	open (O, '>'.$output_file) || die;
	print O "SNP\tAllele\tChromosome\tPos\tCategory\n";	
	
	open (In, $input_file) || die;
	my @snp_pos;
	my %hash_snp;
	my %hash;
#	my $k = 0;
	while (<In>) {
		chomp;
		next if /SNP/;
#		next if $k > 1000;
#		$k++;
		my @t = split /\t/;
		push @snp_pos, $t[3];
		$hash_snp{$t[3]} = $t[0]."\t".$t[1]."\t".$t[2]."\t".$t[3];
		$hash{$t[3]} = $t[4];
	}
	close In;
###########################Use two foreach loops to add the hash tag###################

	foreach my $gene_pos (@$gene_pos_array1) {
		my ($w1_s, $w1_e) = ($gene_pos-$w1, $gene_pos);
		my ($w2_s, $w2_e) = ($gene_pos-$w2, $gene_pos-$w1);
		foreach my $snp_p (@snp_pos) {
			next if ($hash{$snp_p} eq 'P1' || $hash{$snp_p} eq 'P5');
			if ($snp_p >= $w1_s && $snp_p <= $w1_e ){
				$hash{$snp_p} = 'P1';
				} elsif( $snp_p >= $w2_s &&$snp_p < $w2_e){
					$hash{$snp_p} = 'P5';
				} else {
					$hash{$snp_p} = 'Non';
				}	
		}
	}
	
	foreach my $gene_p (@$gene_pos_array2) {
		my ($w1_s, $w1_e) = ( $gene_p, $gene_p+$w1);
		my ($w2_s, $w2_e) = ($gene_p+$w1, $gene_p+$w2);
		foreach my $snp_p (@snp_pos) {
			next unless $hash{$snp_p} eq 'Non';
			if ($snp_p >= $w1_s && $snp_p <= $w1_e ){
				$hash{$snp_p} = 'P1';
				} elsif( $snp_p > $w2_s &&$snp_p <= $w2_e){
					$hash{$snp_p} = 'P5';
				} else {
					$hash{$snp_p} = 'Non';
				}	
		}
	}
	
####################################Print out the result################################################
	foreach my $snp_p (@snp_pos){
		print O $hash_snp{$snp_p}."\t".$hash{$snp_p}."\n";
		}
}

sub Get_gene_pos {
	my ($file, $chro) = @_;
	open (F, $file) || die;
	my (@gene1, @gene2, @gene1_sorted, @gene2_sorted );
	while (<F>) { 
		chomp;
		next if /#/;
		my @t = split /\t/;
		next unless ($t[0] eq $chro && $t[2] eq 'gene'); 
		my $strand = $t[6];
		if ($strand =~ /\+/){
			push @gene1, $t[3];			
			}
		else {
			push @gene2, $t[4];	
			}
		}
	close F;
	@gene1_sorted = sort {$a <=> $b} @gene1;
	@gene2_sorted = sort {$a <=> $b} @gene2;
	return (\@gene1_sorted, \@gene2_sorted);
	}