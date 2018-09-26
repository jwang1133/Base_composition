use strict;
use Bio::DB::Fasta;


#####for later on convenience, add those several parameters####
my $group = 'Maize';
my $group_number=100;
my $missing_percent = 20;
my $maf_percent = 5;
my $group_info = $group.$group_number;
my ($ch_s, $ch_e) = (1, 10);
my ($r_s, $r_e) = ($ARGV[0], $ARGV[1]);
#############Specifying direcotories############
my $dir = '/XXX/';
my $ref_dir = $dir.$group_info.'/XXX/';

my $input_file_dir = $dir.$group_info.'/geno/';
my $output_file_dir = $dir.$group_info.'/trp_motif/';
my $at_chro_size_hashref = Parse_chro_length($dir.$group_info.'/Genome_db/'.$group_info.'_chromosome_length');


my $segment_size = 1000;
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
	
my $output_file = $output_file_dir.$group_info.'_genomewise_random_trp_motif'.'_maf'.$maf_percent.'_miss'.$missing_percent.'r_'.$r_s.'_'.$r_e; 
open (O, '>'.$output_file) || die;
print O "Round\tMotif\tCount\tFrequency\n";


######################calculate the motifs with a loop
my %hash_motif;
my %hash_SNP;
for (my $r = $r_s; $r <= $r_e;  $r ++ ){
	for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	
		my $ch_fas_file = $ref_dir.'Chr'.$ch.'.fa';
		next unless -e $ch_fas_file;
		my $ch_fas_db = Bio::DB::Fasta->new($ch_fas_file);
		my $chro = 'Chr'.$ch;	
		my $input_file = $input_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 
		next unless -e $input_file;

		my $chro_size = $$at_chro_size_hashref{'Chr'.$ch};
		open (Hmp, $input_file) || die;
		my $site_0 = 0;
	
		while (<Hmp>){
		chomp;
		my @t = split /\t/;
		next if $t[1] eq 'alleles';
		my $site = $t[3];
		my $flag = $site - $site_0;
		next if $flag <= 3; 
		   $site_0 = $site;
		my ($a1, $a2) = split /\//, $t[1];
		my $a1a2 = join '', sort($a1, $a2);
		my ($site_p, $site_a) = ($site - 1, $site + 1);
		my $ref_trp_fas = $ch_fas_db -> seq($ch, $site_p => $site_a);
		next if $ref_trp_fas  =~ /N/;
		my $site_t;
		if($site <= $segment_size + 2){$site_t = $site + int(rand(1000));}
		elsif($site > ($segment_size + 2) && $site <= ($chro_size - $segment_size -2)){
			my $random_num = rand();
			if ($random_num <= 0.5){$site_t = $site + int(rand(1000)) ;}
			else {$site_t = $site - int(rand(1000));}
			}
		else {$site_t = $site - int(rand(1000));}
		my ($site_pt, $site_at) = ($site_t - 1, $site_t + 1); 
		my $ref_trp_fas_t = $ch_fas_db -> seq($ch, $site_pt => $site_at);
#		next if $ref_trp_fas  =~ /N/;
		my $ref_site_pt_bp = substr($ref_trp_fas_t,0,1);
		next if $ref_site_pt_bp  =~ /N/;
		my $ref_site_at_bp = substr($ref_trp_fas_t,2,1);
		next if $ref_site_at_bp  =~ /N/;		
		my $di_bp = join('N', $ref_site_pt_bp, $ref_site_at_bp );		
		my $motif_type = join ('_', $a1a2, $di_bp);
		$hash_motif{$r}{$motif_type} ++;
		$hash_SNP{$r} ++;
		}
	close Hmp;
	}

##############print out the motif frequency to the output file
foreach my $mot_t (@motif_type) {
		print O $r."\t".$mot_t."\t";
		my $cnt = exists $hash_motif{$r}{$mot_t} ? $hash_motif{$r}{$mot_t} : 0;
		print O $cnt."\t"; 
		my $total_count = exists $hash_SNP{$r} ? $hash_SNP{$r} : 0;
		my $cnt_freq = sprintf "%.5f", $cnt / $total_count;
		print O $cnt_freq."\n";	
	}
}


sub Parse_chro_length {
	my ($file) = @_;
	my %hash;
	open (FILE, $file) || die;
	while (<FILE>) {
		chomp;
		my @t = split /\t/;
		   $hash{$t[0]} = $t[1];
		}
		close FILE;
		return \%hash;
	}	