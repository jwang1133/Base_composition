use strict;
use Bio::DB::Fasta;

my $ref_dir ='/XXX/';
my $ori_dir = '/XXX/';

my $Mt_db = '/XXX/';

my $out_put_file_dir1 = '/XXX/';
mkdir $out_put_file_dir1 unless -e $out_put_file_dir1;

my $out_put_file_dir2 = '/XXX/';
mkdir $out_put_file_dir2 unless -e $out_put_file_dir2;

my $parameter = "\"m D\"";
my $segment_size = 2000000;

##################################Part that need to be changed, sot hat we can run parallelly 

my ($ch_s, $ch_e, $seg_s, $seg_e) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

#########################################################################################

my $at_chro_size_hashref = Parse_chro_length();

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $ch_fas_file = $ref_dir.'Gm'.$ch.'.fa';
	my $ch_fas_db = Bio::DB::Fasta->new($ch_fas_file);
	my $chro_size = $$at_chro_size_hashref{'Chr'.$ch};	
	
	for (my $seg = $seg_s; $seg <= $seg_e; $seg = $seg+2 ){

		my $ori_seg_file = $ori_dir.'Gm302_chr'.$ch.'_'.$seg.'Mb';
		next unless -e $ori_seg_file;
		my $blast_result = $out_put_file_dir1.'Gm302_chr'.$ch.'_'.$seg.'Mb.snp_blast_result';
		my $seg_ancestral_allele_file = $out_put_file_dir2.'Gm302_Chr'.$ch.'_'.$seg.'Mb.ancestral_allele';
		open (F, $ori_seg_file) || die;
		open (OUTPUT, '>'.$seg_ancestral_allele_file) || die;
		while (<F>) {
			chomp;
			my @t = split /\t/;		
			next if $t[0] =~ /Chromosome/;	
			my $site = $t[1];	
			my $output_file = $out_put_file_dir1.'Gm302_chr'.$ch.'_'.$seg.'Mb.flanking_seq';
			my ($site_p, $site_a) = ($site - 29, $site + 29);	
			next unless (($site_p >= 0) && ($site_a <= $chro_size));
			my $ref_seq_fas = $ch_fas_db -> seq($ch, $site_p => $site_a);
			my $snp_flanking_region_file = Parse_flanking_seq_file($ch, $site, $ref_seq_fas, $output_file);
			system("blastall -p blastn -d  $Mt_db -i $snp_flanking_region_file -m 4 -e 0.1 -o $blast_result  -F $parameter -b 10");
			open (IN, $blast_result) || die;	  
			my $line_num = 0;
			my $flag = 1;
			my $pline_n = 1;
			my $non_base_len;
			my $all_count_len;
			my $snp_position;
			my $ancestral_allele;
			my $flag_b = 0;  
			while(<IN>){
				chomp;
				my $line = $_; 
				$line_num ++;
				next unless ($line_num >=17);
				if ($line =~ /No hits found/){
					$flag = 0;	
				};		
				next unless $flag != 0;	
				
				if ($line =~ /1_0/){
					$pline_n = $line_num;
					$flag ++;
					my $count = 0;
					my @t1 = split //, $line;
					my @t2 = split /\s+/, $line;
					##This is to make sure the matched region including the SNPs
					if (($t2[1] <= 30) && ($t2[3] >= 30)){
						$flag_b = 1;
						my $flag_a = 1;
						##This is to count the pos from the start to the sequence start
						foreach my $e (@t1){
							my $eb = $e; 				
							if($eb =~ /[acgtn]/){
								$flag_a = 0;
							}
							next unless $flag_a != 0;
							if ($eb  !~ /[acgtn]/){
								$count++;
							}
						}
						$non_base_len = $count;
						my @qseq = split//, $t2[2];
						my $all_count = 0;
						my $b_count = 0;
						foreach my $b (@qseq){
							next unless $b_count < 30 - $t2[1] + 1;
							if ($b =~ /[acgtn]/){
								$b_count ++;
							}
							$all_count ++;
						}
						$all_count_len = $all_count;
						$snp_position = $non_base_len + $all_count_len;	 
#						print $snp_position."\n"; 			
					}
				}	
				next unless (($flag == 2)&& ($flag_b == 1)); ##to control incase give multiple sets alignment.	 	 	
			 	
				if ($line_num - $pline_n == 1){
					my @t3 = split /\s+/, $line;
					my $coverage = abs($t3[3] - $t3[1]) + 1;
					$ancestral_allele = substr($line,$snp_position-1,1);
					next unless $ancestral_allele =~ /[acgt]/;
					$ancestral_allele = uc $ancestral_allele;
					print OUTPUT 	$ch."\t".$site."\t".$ancestral_allele."\t".$coverage."\n";
				}	
			}
			close IN;
		}
		close F;
 }
}



sub Parse_chro_length {
	my $file = '/hdd/jinyuw/BCS/Gm302/Genome_db/Gm302_chromosome_length_base_fa'; 
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


sub Parse_flanking_seq_file{
#	my $out_put_file = '/media/lixr/_media_disk_1_/Jinyu/Gm302/snp_flanking_seq/snp_flanking_seq';
	my ($chromosome, $bp, $flanking_sequence, $out_put_file) = @_;
	open (OUT, '>'.$out_put_file) || die;
	print OUT ">chr".$chromosome.'snp'.$bp."\n";
 	print OUT $flanking_sequence."\n";
 	close OUT;
 	return $out_put_file;
 	
#################################################################################################
##there is something wrong with this part, and rememeber that we could not return a filehandle, 
##we can only return a file name that contains the content we want.
#	open my $fh, '>', $out_put_file, or die $!;
#  print $fh ">chr".$chromosome.'snp'.$bp."\n";   
#  print $fh $flanking_sequence."\n";
#  close $fh;
#  return $fh;
###################################################################################################
	}
	
	
	
