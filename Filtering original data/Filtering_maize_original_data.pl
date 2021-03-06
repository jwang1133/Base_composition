#!/usr/bin/perl -w
use strict;

##Specify the input and output file directory
my $input_file_dir = '/XXX/';
my $output_file_dir = '/XXX/';
mkdir $output_file_dir unless -e $output_file_dir; 
my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

##Spefify the minor allele frequency threshold and the missing rate threshold
my $maf_r = 0.05;
my $missing_r = 20;  
my $group_num = 100;

my $missing_rate = $missing_r*100;
my $maf = $maf_r*100;

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $input_file = $input_file_dir.'maizeHapMapV2_B73RefGenV2_201203028_chr'.$ch.'.hmp.txt';
	next unless -e $input_file;
	my $output_file = $output_file_dir.'maize'.$group_num.'_Chr'.$ch.'.maf'.$maf.'_miss'.$missing_rate;
	next if -e $output_file;
	my ($bp_arrayref, $bp_snp_info_hashref) = Parse_chr_snp_files($input_file, $maf_r, $missing_r, $group_num);	
	open (OUT, '>'.$output_file) || die;	
	my @filtered_bp = uniq(@$bp_arrayref);
	foreach my $b (@filtered_bp) {
		next unless exists $$bp_snp_info_hashref{$b};
		my $snp_info_line = $$bp_snp_info_hashref{$b};
		print OUT $snp_info_line."\n";
		}
	}

sub Parse_chr_snp_files{
		my ($input_f, $maf_rate, $missing_f, $group_number) = @_;
		my %bp_snp_info;
		my @bp_array;
		open (F, $input_f) || die;
		while (<F>) {
			chomp;	
			my @t = split /\t/;	
			next if /#/;
			my $alt_ref_allele = $t[1];
			next if $alt_ref_allele =~ /\+/;
			my $bp = $t[3];
			my @snp_info = @t[0..8];
			my @strain_allele;
			for ( my $index =11;  $index <=$#t; $index=$index+1){
			my $i = $index;
			#####remove four accessions that have smaller number of snps
			next if $i == 62;
			next if $i == 85;
			next if $i == 87;
			next if $i == 89;
			##########################################################
			my $genotype = $t[$index];
			my $gg;
			if ($genotype =~ /[ACGTN]/){
			$gg = $genotype;					
				}
			else {$gg = 'N';}
			push @strain_allele, $gg;
			}
			my $total=0;  
			my %hash;
			foreach my $a (@strain_allele) {
				next if $a =~ /N/;
				$total ++;
				$hash{$a} ++;
			}
			my $missing_rate = sprintf "%.2f", ($group_number - $total)/$group_number;
			my @sorted_allele	= sort { $hash{$a} <=> $hash{$b} } keys %hash;
			my $maf = sprintf "%.2f", $hash{$sorted_allele[0]} / $total;
			next unless ($#sorted_allele == 1 && $maf >= $maf_rate && $missing_rate <= $missing_f);
			push @bp_array, $bp;
			push @snp_info, $missing_rate;
			push @snp_info, $maf;
			push @snp_info, @strain_allele;
			my $snp_string = join("\t", @snp_info);
			next if exists $bp_snp_info{$bp};   
			$bp_snp_info{$bp} = $snp_string;
		}
		close F;	
		return (\@bp_array, \%bp_snp_info);	
}




sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}