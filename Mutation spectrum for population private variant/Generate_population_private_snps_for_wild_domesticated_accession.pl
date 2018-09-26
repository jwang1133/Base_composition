use strict;

my $dir = '/XXX/';
my $group_info = 'Gm302';

####################################################
my $Wild_group = 'GS';
my $Wild_group_num = 62; 


my $Dom_group = 'GD';
my $Dom_group_num = 240;

####################################################
my $Wild_group_info = $Wild_group.$Wild_group_num;
my $Dom_group_info = $Dom_group.$Dom_group_num;

my ($ch_s, $ch_e) = (1, 20);

my $input_file_dir = $dir.$group_info.'/XXX/';
my $output_file_dir = $dir.$group_info.'/XXX/';
mkdir $output_file_dir unless -e $output_file_dir;
my $ancestral_allele_dir = $dir.$group_info.'/XXX/';

my ($Wild_index_arrayref,$Wild_id_arrayref) = Get_strain_index_id1($dir.$group_info.'/pop_specific_geno/'.$Wild_group_info.'_from_Gm302_Index');
my ($Dom_index_arrayref,$Dom_id_arrayref) = Get_strain_index_id2($dir.$group_info.'/pop_specific_geno/'.$Dom_group_info.'_from_Gm302_Index');

my @hmp_pre = qw(rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID);

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	
	my $input_file = $input_file_dir.'Gm'.$ch.'.snp.genotype_v2_sorted';
	next unless -e $input_file;
	######################Different ancestral allele genotype file########
	my $ancestral_allele_file = $ancestral_allele_dir.'Gm302_Chr'.$ch.'.ancestral_allele';
	next unless -e $ancestral_allele_file;
	my $ancestral_allele_ref = Parse_ancestral_allele_info($ancestral_allele_file);

	my $Wild_private_file = $output_file_dir. $Wild_group_info.'_specific_Chr'.$ch.'.hmp'; 
	my $Dom_private_file = $output_file_dir.$Dom_group_info.'_specific_Chr'.$ch.'.hmp';
	
	open (OUT1, '>'.$Wild_private_file)|| die;		
	print OUT1 $_."\t" foreach (@hmp_pre);
	print OUT1 'QCcode';
	print OUT1 "\t".$_ foreach (@$Wild_id_arrayref);
	print OUT1 "\n";
	
	
	open (OUT2, '>'.$Dom_private_file)|| die;
	print OUT2 $_."\t" foreach (@hmp_pre);
	print OUT2 'QCcode';
	print OUT2 "\t".$_ foreach (@$Dom_id_arrayref);
	print OUT2 "\n";
 
 
 
 	open (F, $input_file) || die;
	while (<F>) {
			chomp;
			my @t = split /\t/;		
			next if $t[0] =~ /alleles/;	
			my $bp = $t[1];
			next unless exists $$ancestral_allele_ref{$bp};
#			my $ancestral_allele = $$ancestral_allele_ref{$bp};
			my $ancestral_allele_infor = $$ancestral_allele_ref{$bp};
			my @ancestral_array = split(/\t/, $ancestral_allele_infor);
			my ($ancestral_allele, $ancestral_coverage, $ancestral_evalue)  = ($ancestral_array[0], $ancestral_array[1],  $ancestral_array[2]);
			next unless $ancestral_allele  =~ /[ACGT]/;
#			my @allele_info2 = ("+", "RefGenV2", "Gm302", "NA");
			my @allele_info2 = ("+", "RefGenV2", "Gm302");
			push @allele_info2, $ancestral_evalue;
			push @allele_info2, $ancestral_coverage;
			my @allele_info;
			my $pos_x = sprintf "%08d", $bp;
			my $rs_id = 'rs'.$ch.$pos_x; 
			push @allele_info, $rs_id;
			
			my @strain_allele;
			for ( my $index =2;  $index <=$#t; $index=$index+1){
				my $genotype = $t[$index];
				my $gg;

			if ($genotype =~ /[ACGTN]/){
				$gg = $genotype;					
				}
			else {$gg = 'N';}
			push @strain_allele, $gg;
			}																																						  #use this first loop to modify the snp call in the original snp_array
			my $total=0;  
			my %hash;
			foreach my $a (@strain_allele) {
				next if $a =~ /N/;
				$total ++;
				$hash{$a} ++;
			}			
			my @sorted_allele	= sort { $hash{$a} <=> $hash{$b} } keys %hash;
			next unless $#sorted_allele == 1;																								#make sure in the original snp_array, it is a biallelic snp
			my $alt_allele;																						
			if ($sorted_allele[1] eq $ancestral_allele )	{ $alt_allele = $sorted_allele[0]; }
			else {$alt_allele = $sorted_allele[1];}
			my $ancestral_alt_allele = join('/', $ancestral_allele, $alt_allele);		
			push @allele_info, $ancestral_alt_allele;
			push @allele_info, $ch;
			push @allele_info, $bp;
			push @allele_info, @allele_info2;

##################################################################################The above are for get the overal line info
#############Processing the Wild file part######################################################################################			
			my @Wild_bp_allele;																														
			foreach my $Wild_i (@$Wild_index_arrayref)	{
				my $Wild_geno = $strain_allele[$Wild_i-1];
				push @Wild_bp_allele, $Wild_geno;	
			}				
			my $Wild_count_ancestal=0;  
			my $Wild_count_nonancestral=0;				
			my $Wild_total_count = 0;
			my %Wild_hash;
			foreach my $Wild_a (@Wild_bp_allele) {
				next if $Wild_a =~ /N/;
				$Wild_total_count ++;
				$Wild_hash{$Wild_a} ++;
				if ($Wild_a eq $ancestral_allele ){$Wild_count_ancestal++;}
				else {$Wild_count_nonancestral++;}			
			}	
			
		  my $Wild_missing_rate = sprintf "%.2f", ($Wild_group_num - $Wild_total_count)/$Wild_group_num;
			my @Wild_sorted_allele	= sort { $Wild_hash{$a} <=> $Wild_hash{$b} } keys %Wild_hash;
			my $Wild_maf = 0;
			if ($Wild_total_count > 0){
			$Wild_maf = sprintf "%.2f", $Wild_hash{$Wild_sorted_allele[0]} / $Wild_total_count;
		  }
		  
#############Processing the Dom file part######################################################################################			
		  
			my @Dom_bp_allele;
			foreach my $Dom_i (@$Dom_index_arrayref)	{
				my $Dom_geno = $strain_allele[$Dom_i-1];
				push @Dom_bp_allele, $Dom_geno;	
				}	
			my $Dom_count_ancestal=0;  
			my $Dom_count_nonancestral=0;			
			my $Dom_total_count = 0;
			my %Dom_hash;
			foreach my $Dom_a (@Dom_bp_allele) {
				next if $Dom_a =~ /N/;
				$Dom_total_count ++;
				$Dom_hash{$Dom_a} ++;
				if ($Dom_a eq $ancestral_allele ){$Dom_count_ancestal++;}
				else {$Dom_count_nonancestral++;}			
			}							
		  my $Dom_missing_rate = sprintf "%.2f", ($Dom_group_num - $Dom_total_count)/$Dom_group_num;
			my @Dom_sorted_allele	= sort { $Dom_hash{$a} <=> $Dom_hash{$b} } keys %Dom_hash;
			my $Dom_maf = 0;
			if ($Dom_total_count >0){
			$Dom_maf = sprintf "%.2f", $Dom_hash{$Dom_sorted_allele[0]} / $Dom_total_count;		
			}



#############Print out  file part######################################################################################			

			if( $Wild_count_nonancestral != 0 && $Dom_count_nonancestral == 0 ){	
										
			print OUT1 $_."\t" foreach (@allele_info);
			print OUT1 $Wild_missing_rate."\t".$Wild_maf;
			print OUT1 "\t".$_ foreach (@Wild_bp_allele);	
			print OUT1 "\n";
			
			}

			if($Wild_count_nonancestral == 0 && $Dom_count_nonancestral != 0 ){	
										
			print OUT2 $_."\t" foreach (@allele_info);
			print OUT2 $Dom_missing_rate."\t".$Dom_maf;
			print OUT2 "\t".$_ foreach (@Dom_bp_allele);	
			print OUT2 "\n";
			
			}
			
		}	
		close F;
			
	}





	sub  Get_strain_index_id1{
	my $file1 = shift;
	my (@Wild_index, @Wild_strains);
	open (F1, $file1 ) || die;
	while (<F1>) {
		chomp;
		next if /index/;
		my @t = split /\t/;
		push @Wild_index, $t[0];
		push @Wild_strains, $t[1];
		}
	close F1;
	####processing second files
	return (\@Wild_index, \@Wild_strains);
	
}


	sub  Get_strain_index_id2{
	my $file2 = shift;
	my (@Dom_index, @Dom_strains);
	####processing second files
	open (F2, $file2 ) || die;
	while (<F2>) {
		chomp;
		next if /index/;
		my @t = split /\t/;
		push @Dom_index, $t[0];
		push @Dom_strains, $t[1];
		}
	close F2;	
	return (\@Dom_index, \@Dom_strains);	
}


	sub Parse_ancestral_allele_info{
	my $file = shift;
	my %hash_ancestral;
	open (F, $file ) || die;
	while (<F>) {
		chomp;
		next if /Chromosome/;
		my @t = split /\t/;
		$hash_ancestral{$t[1]} = $t[2]."\t".$t[3]."\t".$t[4];
		}
	close F;	
	return \%hash_ancestral;	
	}