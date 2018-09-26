use strict;
use Bio::DB::Fasta;

#################specifying parameters
my $group = 'Maize';
my $group_number=100;
my $missing_percent = 20;
my $maf_percent = 5;
my $group_info = $group.$group_number;
my ($ch_s, $ch_e) = (1, 10);
#my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);

########################specifying directory
my $dir = '/XXX/';
my $ref_dir = $dir.$group_info.'/XXX/';
my $input_file_dir = $dir.$group_info.'/geno/';
my $methyl_file_dir = $dir.$group_info.'/Methylation_analysis/B73RefV4/';
my $genic_non_genic_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';

my $output_file_dir = $dir.$group_info.'/Methylation_analysis/trp_motif/';
mkdir $output_file_dir unless -e $output_file_dir; 		

my $cent_infor_hashref = Parse_cent_infor($dir.$group_info.'/Genome_db/'.$group_info.'_pericent_cent_infor');
#my $peri_segment = 40000000;	


my @bases = qw /A C G T/;
#my @subs_type = qw/AC AG AT CG CT GT TG TC TA GC GA CA/;
my @subs_type = qw/AC AG AT CG CT GT/;
my @di_fas;
foreach my $b1 (@bases) {
	foreach my $b2 (@bases) {
			my $di = $b1.'N'.$b2;
			push @di_fas, $di;		
		}
	}
	
my @motif_type;
foreach my $sub(@subs_type){
	foreach my $di_t(@di_fas){
		my $motif_t =$sub.'_'.$di_t;
		push @motif_type, $motif_t;
		}
	}

my @category_array	= ("peri", "nonperi", "genic", "nongenic", "methyl", "unmethyl", "methyl_peri", "methyl_nonperi", "unmethyl_peri", "unmethyl_nonperi", "methyl_genic", "methyl_nongenic", "unmethyl_genic", "unmethyl_nongenic");
my @motif_type1 = qw/ AG_CNA  AG_CNG  CT_CNG CT_TNG/;
my %hash;
my %hash_cnt;

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $ch_fas_file = $ref_dir.'Chr'.$ch.'.fa';
	next unless -e $ch_fas_file;
	my $ch_fas_db = Bio::DB::Fasta->new($ch_fas_file);
	
	my $cent_infor = $$cent_infor_hashref{$ch};
	my ($pericent_s, $pericent_e) = split /\t/, $cent_infor;
	
	my $input_file = $input_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 	
	next unless -e $input_file;
#	
	my $nongenic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic';
	next unless -e $nongenic_snp_file;
	my $nongenic_snp_info_hashref = Get_nongenic_snp_infor($nongenic_snp_file);

	my $methyl_snp_file = $methyl_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent.'_methyl_infor_sorted';
	next unless -e $methyl_snp_file;
	my $methyl_snp_info_hashref = Get_methyl_snp_infor($methyl_snp_file);
	
	open (Hmp, $input_file) || die;
	my $site_0 = 0;
#  my $k = 0;
	while (<Hmp>){
		chomp;
		my @t = split /\t/;
		next if /alleles/;
		my $site = $t[3];	
		###########################
#		$k ++;
#		last if $k > 10000;
		###############################	
		my $flag = $site - $site_0;
		next if $flag <= 3; 
		$site_0 = $site;    
		my ($a1, $a2) = split /\//, $t[1];
		my $a1a2 = join '', sort($a1, $a2);
		my ($site_p, $site_a) = ($site - 1, $site + 1);
		my $ref_trp_fas = $ch_fas_db -> seq($ch, $site_p => $site_a);
#		next if $ref_trp_fas  =~ /N/;
		my $ref_site_p_bp = substr($ref_trp_fas,0,1);
		next if $ref_site_p_bp  =~ /N/;
		my $ref_site_a_bp = substr($ref_trp_fas,2,1);
		next if $ref_site_a_bp  =~ /N/;
		my $di_bp = join('N', $ref_site_p_bp, $ref_site_a_bp );	
		########Check which category it belongs to
		my ($flag1, $flag2, $flag3, $flag4, $flag5, $flag6, $flag7, $flag8, $flag9, $flag10, $flag11, $flag12, , $flag13, $flag14 ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
		
	  if ($site >= $pericent_s && $site <= $pericent_e){ $flag1 = 1;}
	  else { $flag2 = 1; }
	  
	  if (exists $$nongenic_snp_info_hashref{$site} ){ $flag4 = 1; }
	  else { $flag3 = 1; }
	 
	  if ( exists $$methyl_snp_info_hashref{$site}) {$flag5 = 1; }
	  else {$flag6 = 1; }
	  
	  if ($flag5 == 1 && $flag1 == 1) { $flag7 = 1 ;}
	  if ($flag5 == 1 && $flag2 == 1) { $flag8 = 1 ;}
	  if ($flag6 == 1 && $flag1 == 1) { $flag9 = 1 ;}
	  if ($flag6 == 1 && $flag2 == 1) { $flag10 = 1 ;}
	  
	  if ($flag5 == 1 && $flag3 == 1) { $flag11 = 1 ;}
	  if ($flag5 == 1 && $flag4 == 1) { $flag12 = 1 ;}
	  if ($flag6 == 1 && $flag3 == 1) { $flag13 = 1 ;}
	  if ($flag6 == 1 && $flag4 == 1) { $flag14 = 1 ;}
	  
		for (my $i = 11; $i <= $#t; $i ++) {		
			my $g = $t[$i];
			next if $g =~ /N/;
			my $motif_type = join ('_', $a1a2, $di_bp);
			if ($flag1 == 1) {$hash{$category_array[0]}{$i}{$motif_type} ++; $hash_cnt{$category_array[0]}{$i} ++;}
			if ($flag2 == 1) {$hash{$category_array[1]}{$i}{$motif_type} ++; $hash_cnt{$category_array[1]}{$i} ++;}
			if ($flag3 == 1) {$hash{$category_array[2]}{$i}{$motif_type} ++; $hash_cnt{$category_array[2]}{$i} ++;}
			if ($flag4 == 1) {$hash{$category_array[3]}{$i}{$motif_type} ++; $hash_cnt{$category_array[3]}{$i} ++;}
			if ($flag5 == 1) {$hash{$category_array[4]}{$i}{$motif_type} ++; $hash_cnt{$category_array[4]}{$i} ++;}
			if ($flag6 == 1) {$hash{$category_array[5]}{$i}{$motif_type} ++; $hash_cnt{$category_array[5]}{$i} ++;}
			if ($flag7 == 1) {$hash{$category_array[6]}{$i}{$motif_type} ++; $hash_cnt{$category_array[6]}{$i} ++;}
			if ($flag8 == 1) {$hash{$category_array[7]}{$i}{$motif_type} ++; $hash_cnt{$category_array[7]}{$i} ++;}
			if ($flag9 == 1) {$hash{$category_array[8]}{$i}{$motif_type} ++; $hash_cnt{$category_array[8]}{$i} ++;}
			if ($flag10 == 1) {$hash{$category_array[9]}{$i}{$motif_type} ++; $hash_cnt{$category_array[9]}{$i} ++;}
			if ($flag11 == 1) {$hash{$category_array[10]}{$i}{$motif_type} ++; $hash_cnt{$category_array[10]}{$i} ++;}
			if ($flag12 == 1) {$hash{$category_array[11]}{$i}{$motif_type} ++; $hash_cnt{$category_array[11]}{$i} ++;}
			if ($flag13 == 1) {$hash{$category_array[12]}{$i}{$motif_type} ++; $hash_cnt{$category_array[12]}{$i} ++;}
			if ($flag14 == 1) {$hash{$category_array[13]}{$i}{$motif_type} ++; $hash_cnt{$category_array[13]}{$i} ++;}
		}
	}
	close Hmp;
}

##print out the file	
	my $output_file = $output_file_dir.$group_info.'_genomewise_methyl_unmethyl_peri_nonperi_genic_nongenic_trp_motif'.'_maf'.$maf_percent.'_miss'.$missing_percent; 
	open (O, '>'.$output_file) || die;
	print O "Category\tStrain";
	print O "\t".$_ foreach (@motif_type1);
	print O "\tTotal\n";
			
foreach my $categ (@category_array) {
	for (my $i = 11; $i <= $group_number+10; $i ++) {
		my $ind = $i-10;
		print O $categ."\t".'maize'.$ind;
		foreach my $mot_t (@motif_type1) {
#			print $mot_t."\n";
			my $cnt1 = exists $hash{$categ}{$i}{$mot_t} ? $hash{$categ}{$i}{$mot_t} : 0;
			print O "\t".$cnt1; 
			}
		my $total_cnt = exists $hash_cnt{$categ}{$i}? $hash_cnt{$categ}{$i} : 0;
		print O "\t".$total_cnt;
		print O "\n";	
	}
}	

	
	sub Get_methyl_snp_infor {
	my ($methyl_f) = @_;
	my %snp_infor;
		open (F2, $methyl_f ) || die;
		while (<F2>) {
		chomp;
		next if /Allele/;
		my @t = split /\t/;
		next if $t[4] eq 'unmethyl';
		$snp_infor{$t[3]} = $t[4];  
		}
		close F2;
		return \%snp_infor;
	}
	
	
	sub Get_nongenic_snp_infor {
	my ($nongenic_f) = @_;
	my %snp_infor;
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
	
	
	sub Parse_cent_infor {
	my ($file) = @_;
	my %hash;
	open (FILE, $file) || die;
	while (<FILE>) {
		chomp;
		next if /Chromosome/;
		my @t = split /\t/;
		$hash{$t[0]} = $t[1]."\t".$t[2];
		}
		close FILE;
		return \%hash;
	}	