use strict;
use Bio::DB::Fasta;

#########################specifying parameters
my $group = 'Maize';
my $group_number = 100;
my $missing_percent = 20;
my $maf_percent = 5;
#my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);
my ($ch_s, $ch_e) = (1, 10);

my $group_info = $group.$group_number;

########################specifying directory
my $dir = '/XXXX/';
my $ref_dir = $dir.$group_info.'/XXX/';

my $input_file_dir = $dir.$group_info.'/geno/';
my $output_file_dir = $dir.$group_info.'/trp_motif/';
mkdir $output_file_dir unless -e $output_file_dir;




##########################Define arrays for the 96 type of motifs
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
	


######################calculate the motifs with a loop
my %hash;
	
for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $ch_fas_file = $ref_dir.'Chr'.$ch.'.fa';
	next unless -e $ch_fas_file;
	my $ch_fas_db = Bio::DB::Fasta->new($ch_fas_file);
#	my $chro = 'Chr'.$ch;
	
	my $input_file = $input_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 
	next unless -e $input_file;
	
	open (Hmp, $input_file) || die;
	my $site_0 = 0;

	while (<Hmp>){
		chomp;
		my @t = split /\t/;
		next if $t[1] eq 'alleles';
#		next unless $a1a2 eq 'AG' || $a1a2 eq 'CT';
		my $site = $t[3];
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
		for (my $i = 11; $i <= $#t; $i ++) {		
			my $g = $t[$i];
			next if $g =~ /N/;
			my $motif_type = join ('_', $a1a2, $di_bp);
			$hash{$motif_type} ++;
		}  
	}
	close Hmp;
}	
	
	my $output_file = $output_file_dir.$group_info.'_genomewise_trp_motif'.'_maf'.$maf_percent.'_miss'.$missing_percent; 	
	open (O, '>'.$output_file) || die;
	print O "Motif\tCount\tFrequency";
	print O "\n";		

##Calculate the total count	
	my $total_count = 0;
	foreach my $mot_t (@motif_type) {
			my $cnt = exists $hash{$mot_t} ? $hash{$mot_t} : 0;
			$total_count = $total_count + $cnt;
	}
	
	foreach my $mot_t (@motif_type) {
			print O $mot_t."\t";
			my $cnt = exists $hash{$mot_t} ? $hash{$mot_t} : 0;
			print O $cnt."\t"; 
			my $cnt_freq = sprintf "%.5f", $cnt / $total_count;
			print O $cnt_freq."\n";	
	}
		



