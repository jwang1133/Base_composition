#!/usr/local/bin/perl -w
use Bio::DB::Fasta;
use strict;

####################Specifying the parameters##############
my $group = 'Gm';
my $group_number = 302;
my $missing_percent = 20;
my $maf_percent = 5;
###########################################################
my $missing_rate = $missing_percent/100;
my $maf = $maf_percent/100;
my $group_info = $group.$group_number;

my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);
####################Specifying the directory###############
my $dir = '/XXX/';
my $ori_dir = $dir.$group_info.'/XXX/';
my $hmp_dir = $dir.$group_info.'/XXX/';
mkdir $hmp_dir unless -e $hmp_dir;

my $ref_dir = $dir.$group_info.'/XXX/';
###########################################################

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $hmp_file = $hmp_dir.$group_info.'_Chr'.$ch.'.hmp_maf'.$maf_percent.'_miss'.$missing_percent; 
#	my $chch = $ch < 10 ? '0'.$ch : $ch;
	my $ori_ch_file = $ori_dir.'Gm'.$ch.'.snp.genotype_v2_sorted';
	my $ch_fas_file = $ref_dir.'Gm'.$ch.'.fa';
	my $ch_db =  Bio::DB::Fasta->new($ch_fas_file);
	
	open (Hmp, '>'.$hmp_file) || die;
	open (F, $ori_ch_file) || die;
	
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next if /#/;
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
		my $missing_rate1 = sprintf "%.3f", ($group_number - $total)/$group_number;	
		next unless $missing_rate1 <= $missing_rate; ##This is to control later it will not divide zero of $total
		my @sorted_allele	= sort { $hash{$a} <=> $hash{$b} } keys %hash;
		my $maf_r = sprintf "%.3f", $hash{$sorted_allele[0]} / $total;
		next unless ($#sorted_allele == 1 && $maf_r >= $maf );
		my $major_allele = $sorted_allele[1]; 
		my $minor_allele = $sorted_allele[0]; 
		my $ref_allele = $ch_db->seq($ch, $t[1], $t[1]);
		next unless $ref_allele =~ /[ACGT]/;
		next unless ($ref_allele eq $major_allele || $ref_allele eq $minor_allele);
		my $alt_allele = $ref_allele eq $major_allele ? $minor_allele : $major_allele;
		my $ref_alt_allele = join('/', $ref_allele,$alt_allele);			
		my $pos_x = sprintf "%08d", $bp;
		my $rs_id = 'rs'.$ch.$pos_x;
		print Hmp $rs_id."\t".$ref_alt_allele."\t".$ch."\t".$bp."\t+\tRefGenV2\tGm302\tNA\tNA\t$missing_rate1\t$maf_r";
		print Hmp "\t".$_ foreach (@bp_allele);
		print Hmp "\n";		
			
	}
	close F;	
}
