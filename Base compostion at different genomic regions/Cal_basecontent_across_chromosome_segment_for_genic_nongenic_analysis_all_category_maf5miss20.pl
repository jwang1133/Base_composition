#!/usr/local/bin/perl -w
###This file are modified to calculate the BCS across chromosome segement for 8 subtype (all substype included)
###And only consider three category (Genic, nongenic and all)
use strict;
## merge segment hmps into chr
##those are the parameter need to change


###############Specifying the parameters
my $group = 'Maize';
my $group_number= 100;
my $missing_percent = 20;
my $maf_percent = 5;
###################################
my $group_info = $group.$group_number;

###############Specifying the segment and step infomration
my ($ch_s, $ch_e) = ($ARGV[0], $ARGV[1]);
my $segment_size = 5000000;
my $step_size = 1000000;
my $w = $segment_size - $step_size;
###########################################

my $seg_s = int($segment_size/1000000);
my $step_s = int($step_size/1000000);
my $suffix = $seg_s.'Mb_window_'.$step_s.'Mb_step';

my $dir = '/hdd/jinyuw/BCS/';

my $input_snp_file_dir = $dir.$group_info.'/geno/';
my $genic_non_genic_file_dir = $dir.$group_info.'/Genic_nongenic_analy/';

my $id_line_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_1st_id_line';
my ($strain_id_arrayref, $first_id_line) = Get_strain_id($id_line_file); 

my $strain_group_infor_file = $dir.$group_info.'/hmp_origin/'.$group_info.'_group_infor';
my $strain_group_hashref = Get_group_infor($strain_group_infor_file);

my $chromosome_size_file = $dir.$group_info.'/Genome_db/'.$group_info.'_chromosome_length';
my $at_chro_size_hashref = Parse_chro_length($chromosome_size_file);

#########################################################
#my @subs_type = qw/AC AG AT CG CT GT AG_CT Sub6/;
my @subs_type = qw/Sub6 AG_CT/;
##Sub6 means include all 6 subtype
#my $sub = 'AC';
#my $sub = 'AG';
#my $sub = 'AT';
#my $sub = 'CG';
#my $sub = 'CT';
#my $sub = 'GT';
#########################################################

my @category_array	= ("All", "Genic", "Nongenic","Nonsyn","UTRs","Og","Promoter5k","In", "Syn", "Nyn", "3U", "5U", "Sg", "Rg","Non",  "P5", "P1");

my @category_array1	= ("All", "Genic", "Nongenic","Nonsyn","UTRs","Og","Promoter5k","In", "Syn", "Nyn", "Non" );

my @Genic_array = ("In", "Syn", "Nyn", "3U", "5U", "Sg", "Rg" );
my @Nongenic_array = ("Non", "P5", "P1"); 
####Different ways to classify the Nonsynonymous SNPs: "Nonsyn"#####
my @Nsy_array = ("Nyn", "Sg", "Rg" );
#############################"UTRs"#####################
my @UTRs_array = ( "3U", "5U");
#############################"Og"#####################
my @Othergenic_array = ("Sg", "Rg" );

my @Promoter5K_array = ("P5", "P1");




my @bases = qw /A C G T/;
my (%strain_base_count, %strain_snp_count, %category_subtype_count, %category_count);

my $output_dir = $genic_non_genic_file_dir.'subgroup_analy/';
mkdir $output_dir unless -e $output_dir; 			
my $output_file = $output_dir.$group_info.'_2subgroup_and_across_chro_segment_BCS_all_category_maf'.$maf_percent.'_miss'.$missing_percent.'_across_'.$suffix;
open (OUT, '>'.$output_file) || die;
print OUT "Category\tChromosome\tSegment\tStrain\tGroup\tSubtype\tSub_freq\tA\tC\tG\tT\tTotal\n";
#print OUT "Category\tStrain\tGroup\tSubtype\tSub_freq\tA\tC\tG\tT\tTotal\n";

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	
	my $input_chrsnp_file = $input_snp_file_dir.$group_info.'_Chr'.$ch.'.maf'.$maf_percent.'_miss'.$missing_percent; 	
	next unless -e $input_chrsnp_file;
	my $genic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_genic';
	next unless -e $genic_snp_file;
	my $nongenic_snp_file = $genic_non_genic_file_dir.$group_info.'_Chr'.$ch.'.vcf_maf'.$maf_percent.'_miss'.$missing_percent.'_nongenic_splitted';
	next unless -e $nongenic_snp_file;
	
	my $snp_info_hashref = Get_snp_infor($genic_snp_file, $nongenic_snp_file);
	
	my $chro_size = $$at_chro_size_hashref{'Chr'.$ch};
	my $seg_st = 1;
	my $seg_e = int($chro_size/$w)+1;
#	my $seg_e = 1;
	
		for (my $seg = $seg_st; $seg <= $seg_e ; $seg ++) { 
			my ($segment_s, $segment_e);
			if ($seg == 1){
				$segment_s = 1;
				$segment_e = $segment_size;
				}
			else {
				$segment_s = ($seg-1)*$w;
				$segment_e = ($seg-1)*$w + $segment_size > $chro_size? $chro_size : $seg*$segment_size - $step_size;
				}
			open (INPUT, $input_chrsnp_file ) || die;
			while (<INPUT>) {
				chomp;
				next if /alleles/;
				my @t = split /\t/;
#				$k++;
#				next if $k > 10000;
				my @a1a2 = split /\//,$t[1];
				my $sub = join '', sort @a1a2;
				my $bp = $t[3];
				next unless (($bp>=$segment_s)&&($bp<$segment_e));
					
				my ($flag, $flag1, $flag2, $flag3, $flag4, $flag5, $flag6 ) = qw/0 0 0 0 0 0 0/;
				my $category = $$snp_info_hashref{$bp};
				
				##Calculate the SNP number for that segment for the category
				$category_count{$category_array[0]}{$ch}{$seg}++;
				$category_count{$category}{$ch}{$seg}++;
				
				$category_subtype_count{$category_array[0]}{$ch}{$seg}{$subs_type[0]}++;
				$category_subtype_count{$category}{$ch}{$seg}{$subs_type[0]}++;
				
				##############################################################
				if($sub eq 'AG' || $sub eq 'CT'){
					$flag = 1;
					$category_subtype_count{$category_array[0]}{$ch}{$seg}{$subs_type[1]}++;
					$category_subtype_count{$category}{$ch}{$seg}{$subs_type[1]}++;
				}	
			
#				print $category."\n";
############################################################
				if (grep {$_ eq $category} @Genic_array) {
					$flag1 = 1;
					$category_subtype_count{$category_array[1]}{$ch}{$seg}{$subs_type[0]}++;
					$category_count{$category_array[1]}{$ch}{$seg}++;
					if($flag == 1){
						$category_subtype_count{$category_array[1]}{$ch}{$seg}{$subs_type[1]}++;
					}
				}
##############################4 extra different cases for subtype_frequency######################################3

				if (grep {$_ eq $category} @Nongenic_array) {
					$flag2 = 1;
					$category_subtype_count{$category_array[2]}{$ch}{$seg}{$subs_type[0]}++;
					$category_count{$category_array[2]}{$ch}{$seg}++;
					if($flag == 1){
						$category_subtype_count{$category_array[2]}{$ch}{$seg}{$subs_type[1]}++;
					}
				}
#			
				if (grep {$_ eq $category} @Nsy_array) {
					$flag3 = 1;
					$category_subtype_count{$category_array[3]}{$ch}{$seg}{$subs_type[0]}++;
					$category_count{$category_array[3]}{$ch}{$seg}++;
					if($flag == 1){
						$category_subtype_count{$category_array[3]}{$ch}{$seg}{$subs_type[1]}++;
					}
				}
			
				if (grep {$_ eq $category} @UTRs_array) {
					$flag4 = 1;
					$category_subtype_count{$category_array[4]}{$ch}{$seg}{$subs_type[0]}++;
					$category_count{$category_array[4]}{$ch}{$seg}++;
					if($flag == 1){
						$category_subtype_count{$category_array[4]}{$ch}{$seg}{$subs_type[1]}++;
					}
				}	
				
				if (grep {$_ eq $category} @Othergenic_array) {
					$flag5 = 1;
					$category_subtype_count{$category_array[5]}{$ch}{$seg}{$subs_type[0]}++;
					$category_count{$category_array[5]}{$ch}{$seg}++;
					if($flag == 1){
						$category_subtype_count{$category_array[5]}{$ch}{$seg}{$subs_type[1]}++;
					}
				}	
				
				if (grep {$_ eq $category} @Promoter5K_array) {
					$flag6 = 1;
					$category_subtype_count{$category_array[6]}{$ch}{$seg}{$subs_type[0]}++;
					$category_count{$category_array[6]}{$ch}{$seg}++;
					if($flag == 1){
						$category_subtype_count{$category_array[6]}{$ch}{$seg}{$subs_type[1]}++;
					}
				}	
###################################go through the loop to calculate the base content#################################################								
				for (my $i = 11; $i <= $#t; $i ++) {
					next unless $t[$i] =~ /[ACGT]/;
					my $strain = $$strain_id_arrayref[$i];
					$strain_base_count{$category}{$ch}{$seg}{$sub}{$subs_type[0]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category}{$ch}{$seg}{$sub}{$subs_type[0]}{$strain} ++;		
											
					$strain_base_count{$category_array[0]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
					$strain_snp_count{$category_array[0]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;
							
					if ($flag == 1){
						$strain_base_count{$category}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
													
						$strain_base_count{$category_array[0]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[0]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;
					}

					if ($flag1 == 1) {
						$strain_base_count{$category_array[1]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[1]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;	
						if ($flag == 1){
							$strain_base_count{$category_array[1]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
							$strain_snp_count{$category_array[1]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
							}
					}
					
					if ($flag2 == 1) {
						$strain_base_count{$category_array[2]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[2]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;	
						if ($flag == 1){
							$strain_base_count{$category_array[2]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
							$strain_snp_count{$category_array[2]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
							}
					}
					
					if ($flag3 == 1) {
						$strain_base_count{$category_array[3]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[3]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;	
						if ($flag == 1){
							$strain_base_count{$category_array[3]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
							$strain_snp_count{$category_array[3]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
						}
					}
					
					if ($flag4 == 1) {
						$strain_base_count{$category_array[4]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[4]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;	
						if ($flag == 1){
							$strain_base_count{$category_array[4]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
							$strain_snp_count{$category_array[4]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
						}
					}
					
					if ($flag5 == 1) {
						$strain_base_count{$category_array[5]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[5]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;	
						if ($flag == 1){
							$strain_base_count{$category_array[5]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
							$strain_snp_count{$category_array[5]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
						}
					}
					
					if ($flag6 == 1) {
						$strain_base_count{$category_array[6]}{$ch}{$seg}{$subs_type[0]}{$strain}{$t[$i]} ++;
						$strain_snp_count{$category_array[6]}{$ch}{$seg}{$subs_type[0]}{$strain} ++;	
						if ($flag == 1){
							$strain_base_count{$category_array[6]}{$ch}{$seg}{$subs_type[1]}{$strain}{$t[$i]} ++;
							$strain_snp_count{$category_array[6]}{$ch}{$seg}{$subs_type[1]}{$strain} ++;	
						}
					}
					
				}
				
		}
		close INPUT;	

	}

################################################Print the result out with for loop##########################################
		foreach my $categ (@category_array1) {			
			my @segment = keys %{ $category_subtype_count{$categ}{$ch} };
			foreach my $segment (sort {$a <=> $b} @segment ){
					my $category_total_num = exists $category_count{$categ}{$ch}{$segment} ? $category_count{$categ}{$ch}{$segment} : 1;
					foreach my $sub_t (@subs_type){
						my $category_sub_count = exists $category_subtype_count{$categ}{$ch}{$segment}{$sub_t} ? $category_subtype_count{$categ}{$ch}{$segment}{$sub_t} : 0;
						my $categ_sub_freq = sprintf "%.5f", $category_sub_count/$category_total_num;
#						print OUT $categ."\t".$ch."\t".$segment."\t".$sub_t."\t".$categ_sub_freq."\n";
								for (my $i = 11; $i < @$strain_id_arrayref; $i ++)  {
									my $strain = $$strain_id_arrayref[$i];
									my $group = $$strain_group_hashref{$strain};
									print OUT $categ."\t".$ch."\t".$segment."\t".$strain."\t".$group."\t".$sub_t."\t".$categ_sub_freq;
									my $total_num = exists $strain_snp_count{$categ}{$ch}{$segment}{$sub_t}{$strain} ? $strain_snp_count{$categ}{$ch}{$segment}{$sub_t}{$strain} : 1;
									foreach my $b (@bases) {
										my $base_count = exists $strain_base_count{$categ}{$ch}{$segment}{$sub_t}{$strain}{$b} ? $strain_base_count{$categ}{$ch}{$segment}{$sub_t}{$strain}{$b} : 0;
										my $percent = sprintf "%.5f", $base_count / $total_num;
										print OUT "\t".$percent;
									}
									print OUT "\t".$total_num."\n";
								}	
					}
			}
		}
######################################################################################3
	}
		
		








sub Get_strain_id {
	my $file = shift;
	open (FILE, $file ) || die;
	my (@array, $line);
	while (<FILE>) {
		chomp;
			$line = $_;
		  my @t = split /\t/, $line;
		  @array = @t;  
		}
		close FILE;
		return (\@array, $line);
	}	
	
sub Get_group_infor {
	my $file = shift;
	open (FILE, $file ) || die;
	my %hash;
	while (<FILE>) {
		chomp;
		next if /Group/;
		my @t = split /\t/;
		$hash{$t[0]} = $t[1];  
		}
		close FILE;
		return (\%hash);
	}	

sub Get_snp_infor {
	my ($genic_f, $nongenic_f) = @_;
	my %snp_infor;
		open (F1, $genic_f ) || die;
		while (<F1>) {
		chomp;
		next if /Allele/;
		my @t = split /\t/;
		$snp_infor{$t[3]} = $t[4];  
		}
		close F1;
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
	
	
	
	
	sub Parse_chro_length {
	my $file = shift;
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