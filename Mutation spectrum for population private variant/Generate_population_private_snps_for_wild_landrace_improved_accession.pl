use strict;

my $dir = '/XXX/';

my $group_info = 'Gm302';
####################################################
my $Wild_group = 'GS';
my $Wild_group_num = 62; 

my $Imp_group = 'GC';
my $Imp_group_num = 110;

my $Lan_group = 'GL';
my $Lan_group_num = 130;
####################################################
my ($ch_s, $ch_e) = (1, 20);


my $input_file_dir = $dir.$group_info.'/XXX/';
my $output_file_dir = $dir.$group_info.'/XXX/e1_combined/';
mkdir $output_file_dir unless -e $output_file_dir;
my $ancestral_allele_dir = $dir.$group_info.'/XXX/e1_combined/';



my $Wild_group_info = $Wild_group.$Wild_group_num;
my $Imp_group_info = $Imp_group.$Imp_group_num;
my $Lan_group_info = $Lan_group.$Lan_group_num;


my ($Wild_index_arrayref,$Wild_id_arrayref) = Get_strain_index_id1($dir.$group_info.'/pop_specific_geno/'.$Wild_group_info.'_from_Gm302_Index');
my ($Imp_index_arrayref,$Imp_id_arrayref) = Get_strain_index_id2($dir.$group_info.'/pop_specific_geno/'.$Imp_group_info.'_from_Gm302_Index');
my ($Lan_index_arrayref,$Lan_id_arrayref) = Get_strain_index_id3($dir.$group_info.'/pop_specific_geno/'.$Lan_group_info.'_from_Gm302_Index');

my @hmp_pre = qw(rs# alleles chrom pos strand assembly center protLSID assayLSID panel);

for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	
	my $input_file = $input_file_dir.'Gm'.$ch.'.snp.genotype_v2_sorted';
	next unless -e $input_file;
	######################Different ancestral allele genotype file########
	my $ancestral_allele_file = $ancestral_allele_dir.'Gm302_Chr'.$ch.'.ancestral_allele';
	next unless -e $ancestral_allele_file;
	my $ancestral_allele_ref = Parse_ancestral_allele_info($ancestral_allele_file);

	my $Wild_private_file = $output_file_dir.$Wild_group_info.'_specific_Chr'.$ch.'.hmp'; 
	my $Imp_private_file = $output_file_dir.$Imp_group_info.'_specific_Chr'.$ch.'.hmp';
	my $Lan_private_file = $output_file_dir.$Lan_group_info.'_specific_Chr'.$ch.'.hmp';
	
	open (OUT1, '>'.$Wild_private_file)|| die;		
	print OUT1 $_."\t" foreach (@hmp_pre);
	print OUT1 'QCcode';
	print OUT1 "\t".$_ foreach (@$Wild_id_arrayref);
	print OUT1 "\n";
	
	
	open (OUT2, '>'.$Imp_private_file)|| die;
	print OUT2 $_."\t" foreach (@hmp_pre);
	print OUT2 'QCcode';
	print OUT2 "\t".$_ foreach (@$Imp_id_arrayref);
	print OUT2 "\n";
 
 	open (OUT3, '>'.$Lan_private_file)|| die;
	print OUT3 $_."\t" foreach (@hmp_pre);
	print OUT3 'QCcode';
	print OUT3 "\t".$_ foreach (@$Lan_id_arrayref);
	print OUT3 "\n";
 
 
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

#############Processing the Wild file part######################################################################################			
			my @Wild_bp_allele;																														
			foreach my $Wild_i (@$Wild_index_arrayref)	{
				my $Wild_geno = $strain_allele[$Wild_i - 1];
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
		  
#############Processing the Imp file part######################################################################################			
		  
			my @Imp_bp_allele;
			foreach my $Imp_i (@$Imp_index_arrayref)	{
				my $Imp_geno = $strain_allele[$Imp_i - 1];
				push @Imp_bp_allele, $Imp_geno;	
				}	
			my $Imp_count_ancestal=0;  
			my $Imp_count_nonancestral=0;			
			my $Imp_total_count = 0;
			my %Imp_hash;
			foreach my $Imp_a (@Imp_bp_allele) {
				next if $Imp_a =~ /N/;
				$Imp_total_count ++;
				$Imp_hash{$Imp_a} ++;
				if ($Imp_a eq $ancestral_allele ){$Imp_count_ancestal++;}
				else {$Imp_count_nonancestral++;}			
			}							
		  my $Imp_missing_rate = sprintf "%.2f", ($Imp_group_num - $Imp_total_count)/$Imp_group_num;
			my @Imp_sorted_allele	= sort { $Imp_hash{$a} <=> $Imp_hash{$b} } keys %Imp_hash;
			my $Imp_maf = 0;
			if ($Imp_total_count >0){
			$Imp_maf = sprintf "%.2f", $Imp_hash{$Imp_sorted_allele[0]} / $Imp_total_count;		
			}

#############Processing the Lan file part######################################################################################			
			my @Lan_bp_allele;
			foreach my $Lan_i (@$Lan_index_arrayref)	{
				my $Lan_geno = $strain_allele[$Lan_i - 1];
				push @Lan_bp_allele, $Lan_geno;	
				}	
			my $Lan_count_ancestal=0;  
			my $Lan_count_nonancestral=0;			
			my $Lan_total_count = 0;
			my %Lan_hash;
			foreach my $Lan_a (@Lan_bp_allele) {
				next if $Lan_a =~ /N/;
				$Lan_total_count ++;
				$Lan_hash{$Lan_a} ++;
				if ($Lan_a eq $ancestral_allele ){$Lan_count_ancestal++;}
				else {$Lan_count_nonancestral++;}			
			}							
		  my $Lan_missing_rate = sprintf "%.2f", ($Lan_group_num - $Lan_total_count)/$Lan_group_num;
			my @Lan_sorted_allele	= sort { $Lan_hash{$a} <=> $Lan_hash{$b} } keys %Lan_hash;
			my $Lan_maf = 0;
			if ($Lan_total_count >0){
			$Lan_maf = sprintf "%.2f", $Lan_hash{$Lan_sorted_allele[0]} / $Lan_total_count;		
			}


#############Print out  file part######################################################################################			

			if( $Wild_count_nonancestral != 0 && $Imp_count_nonancestral == 0 && $Lan_count_nonancestral == 0 ){	
										
			print OUT1 $_."\t" foreach (@allele_info);
			print OUT1 $Wild_missing_rate."\t".$Wild_maf;
			print OUT1 "\t".$_ foreach (@Wild_bp_allele);	
			print OUT1 "\n";
			
			}

			if($Wild_count_nonancestral == 0 && $Imp_count_nonancestral != 0 && $Lan_count_nonancestral == 0  ){	
										
			print OUT2 $_."\t" foreach (@allele_info);
			print OUT2 $Imp_missing_rate."\t".$Imp_maf;
			print OUT2 "\t".$_ foreach (@Imp_bp_allele);	
			print OUT2 "\n";
			
			}

			if($Wild_count_nonancestral == 0 && $Imp_count_nonancestral == 0 && $Lan_count_nonancestral != 0  ){	
										
			print OUT3 $_."\t" foreach (@allele_info);
			print OUT3 $Lan_missing_rate."\t".$Lan_maf;
			print OUT3 "\t".$_ foreach (@Lan_bp_allele);	
			print OUT3 "\n";
			
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
	my (@Imp_index, @Imp_strains);
	####processing second files
	open (F2, $file2 ) || die;
	while (<F2>) {
		chomp;
		next if /index/;
		my @t = split /\t/;
		push @Imp_index, $t[0];
		push @Imp_strains, $t[1];
		}
	close F2;	
	return (\@Imp_index, \@Imp_strains);	
}


	sub  Get_strain_index_id3{
	my $file3 = shift;
	my (@Lan_index, @Lan_strains);
	####processing third files
	open (F3, $file3 ) || die;
	while (<F3>) {
		chomp;
		next if /index/;
		my @t = split /\t/;
		push @Lan_index, $t[0];
		push @Lan_strains, $t[1];
		}
	close F3;	
	return (\@Lan_index, \@Lan_strains);	
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