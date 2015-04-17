#This package contains the most frequently used subroutines
#Contact Fan wei, fanw@genomics.org.cn
#Created on 2006-9-20

package  Collect;
use strict qw(subs refs);
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw( Display_seq Complement_Reverse Remask_seq Compare_position Get_Repeats
	Read_fasta Read_psl Read_RepeatMasker);



##blast m8 file format
#Fields: Query id, Subject id, % identity, alignment length, 
#        mismatches, gap openings, q. start, q. end, s. start, s. end, 
#		e-value, bit score


##Definition
##Fields are: <seqname> <source> <feature> <start> <end> <score> <strand> 
##<frame> [attributes] [comments]


##read Genbank format file downlaod from NCBI
#############################################
sub read_Genbank{
	my $file = shift;
	my $hash = shift;
	$/="\nLOCUS";
	open IN, $file || die "fail $file\n";
	while (<IN>) {
		chomp;
		$_ = "LOCUS".$_;
		my ($locus,$len,$type,$class,$ACCESSION,$ORGANISM,$DEFINITION,$ORF,$seq);
		($locus,$len,$type,$class) = ($1,$2,$3,$4) if(/LOCUS\s+(\S+)\s+(\d+)\s+bp\s+(\S+)\s+\S+\s+(\S+)/);
		next unless($locus);
		$ACCESSION = $1 if(/\nACCESSION\s+(\S+)/);	
		$ORGANISM = $1 if(/\n\s+ORGANISM\s+(.+)/);
		$DEFINITION = $1 if(/\nDEFINITION\s+(.+?)\s+ACCESSION/s);
		$DEFINITION =~ s/\n           //g;

		$ORF = $1 if(/\n\s+CDS\s+(\d+\.\.\d+)\s+/); ##complete mRNA sequence
		
		$seq = $1 if(/\nORIGIN(.+?)\/\//s);      
		$seq =~ s/\s+\d+\s+//g;
		$seq =~ s/\s//g;
		
		$hash->{$ACCESSION}{LOCUS} = $locus;
		$hash->{$ACCESSION}{LEN} = $len;
		$hash->{$ACCESSION}{TYPE} = $type;
		$hash->{$ACCESSION}{CLASS} = $class;
		$hash->{$ACCESSION}{ORGANISM} = $ORGANISM;
		$hash->{$ACCESSION}{DEFINITION} = $DEFINITION;
		$hash->{$ACCESSION}{ORF} = $ORF;
		$hash->{$ACCESSION}{SEQ} = $seq;
		$hash->{$ACCESSION}{CONTENT} = $_;
		print $_,"\n";
	}
	close IN;
	$/="\n";
}
##get program dir, derived from LiShengting
#############################################
sub get_bin_path{
	my $BIN_PATH='.';
	my $pwd=`pwd`;
	chomp $pwd;
	
	#print $0,"\n";
	
	if ($0=~/(^\..*)\/.*?/) {
		$BIN_PATH=$pwd."/".$1;
	}elsif ($0=~/(^\/.*)\/.*?/) {
		$BIN_PATH=$1;
	}elsif ($0=~/(.+)\/.*$/) {
		$BIN_PATH=$pwd."/".$1;
	}

	return $BIN_PATH;
}

##static specified character number in a specified string
##usage: my $num = char_num($str,$char);
#############################################
sub char_num{
	my ($str,$char) = @_;
	my $num = 0;
	for (my $i=0; $i<length($str); $i++) {
		$num++ if(substr($str,$i,1) eq $char);
	}
	return $num;
}

#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name}{head} =  $head;
		$hash_p->{$name}{len} = length($seq);
		$hash_p->{$name}{seq} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}


##read psl file
##convert to normal coordinate
#usage: Read_psl($file,\%hash);
####################################################
sub Read_psl{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		s/^\s+//;
		my @temp = split(/\s+/);
		my $qname = $temp[9];
		my $tname = $temp[13];
		my $strand = $temp[8];
		my $gene_start = $temp[15]+1;
		my $gene_end = $temp[16];
		
		my @tstarts=split(/,/,$temp[20]);
		my @tsizes=split(/,/,$temp[18]);

		my (@exon,@intron);
		for (my $i=0; $i<@tstarts; $i++) {
			push @exon, [$tstarts[$i]+1,$tstarts[$i]+$tsizes[$i]];
		}
		for (my $i=0; $i<@tstarts-1; $i++) {
			push @intron, [$tstarts[$i]+$tsizes[$i]+1,$tstarts[$i+1]];
		}

		push @{$hash_p->{$tname}}, [$qname,$strand,$gene_start,$gene_end,\@exon,\@intron];
		
		$total_num++;
	}
	close(IN);

	return $total_num;
}


##read RepeatMasker .out file
#usage: Read_RepeatMasker($file,\%hash);
############################################
sub Read_RepeatMasker{
	my $file=shift;
	my $hash_p=shift; 
	
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		s/^\s+//;
		my @temp=split(/\s+/);
		next if($temp[8] != 'C' && $temp[8] != '+');
		my $tname = $temp[4];
		my $strand = ($temp[8] eq '+') ? '+' : '-';
		my $start = $temp[5];
		my $end = $temp[6];
		my $TE_name = $temp[9];
		my $TE_class = $temp[10];

		push @{$hash_p->{$tname}}, [$TE_name,$strand,$start,$end,$TE_class]; 
	}
	close(IN);

}



#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################

#change the gene in the plus chain to the regular format
#deal with sequence and position automatically
#usage: Complement_Reverse(\$seq);
#		Complement_Reverse(\@pos);
#############################################
sub Complement_Reverse{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
		$$seq_p=~tr/AGCTagct/TCGAtcga/;
		$$seq_p=reverse($$seq_p);  
	}else{ ##deal with position in two demension array
		@$seq_p = reverse @$seq_p;
		foreach my $p (@$seq_p) {
			($p->[0],$p->[1]) = ($p->[1],$p->[0]);
		}

	}
}
#############################################


##remask a sequence with repeat positon 
##usage: Remask_seq(\$seq,\@pos);
###########################################
sub Remask_seq{
	my $seq_p = shift; ##sequence pointer
	my $rep_ap = shift; ##array pointer for repeat position
	my $maskC = (@_) ? shift : "N"; ##set the base for masking result
	
	$$seq_p =~ s/\s//g;
	foreach my $p (@$rep_ap) {
		my ($start,$end) = ($p->[0] <= $p->[1]) ? ($p->[0] , $p->[1]) : ($p->[1], $p->[0]);
		substr($$seq_p,$start-1,$end-$start+1) = $maskC x ($end-$start+1);
	}
}

##get list of repeat starting and ending positons
##process 1 chromosome in a string at a time
##usage: Get_Repeats(\$seq,\@rep);
#############################################
sub Get_Repeats{
	my $seq_p=shift;
	my $rep_ap=shift; #hash pointer
	my $maskC = (@_) ? shift : 'N';  # N or X, that stands for repeat region
	
	$$seq_p =~ s/\s//g;
	while ($$seq_p =~ /$maskC+/g) {
		my $end = pos($$seq_p);
		my $start = $end - length($&) + 1;
		push @$rep_ap, [$start,$end];
	}
}



##conjoin the overlapped fragments, and caculate the redundant size
##usage: conjoin_fragment(\@pos);
##		 my ($all_size,$pure_size,$redunt_size) = conjoin_fragment(\@pos);
##Note than change the pointer's value can cause serious confusion.
sub Conjoin_fragment{
	my $pos_p = shift; ##point to the two dimension input array
	my $new_p = [];         ##point to the two demension result array
	
	my ($all_size, $pure_size, $redunt_size) = (0,0,0); 
	
	return (0,0,0) unless(@$pos_p);

	foreach my $p (@$pos_p) {
			($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
			$all_size += abs($p->[0] - $p->[1]) + 1;
	}
	
	@$pos_p = sort {$a->[0] <=>$b->[0]} @$pos_p;
	push @$new_p, (shift @$pos_p);
	
	foreach my $p (@$pos_p) {
			if ( ($p->[0] - $new_p->[-1][1]) <= 0 ) { # conjoin
					if ($new_p->[-1][1] < $p->[1]) {
							$new_p->[-1][1] = $p->[1]; 
					}
					
			}else{  ## not conjoin
					push @$new_p, $p;
			}
	}
	@$pos_p = @$new_p;

	foreach my $p (@$pos_p) {
			$pure_size += abs($p->[0] - $p->[1]) + 1;
	}
	
	$redunt_size = $all_size - $pure_size;
	return ($all_size,$pure_size,$redunt_size);
}

##compare two position sets, to find overlap region
sub Compare_position{
	my ($pos1_p,$pos2_p) = @_;
	my $pure1_size = Conjoin_fragment($pos1_p)->[1];
	my $pure2_size = Conjoin_fragment($pos2_p)->[1];
	my @pos_comb = (@$pos1_p,@$pos2_p);
	my $share_size = Conjoin_fragment(\@pos_comb)->[2];
	
	print STDERR "pure size of pos1: $pure1_size\n";
	print STDERR "pure size of pos2: $pure2_size\n";
	print STDERR "(1,2) shared size: $share_size\n";
	
	return \@pos_comb;
}

##find overlap region between two fragment sets
##pre-conditon: the two two fragment sets must be no redundant
##position stored in two dimension array.
##usage: conjoin_fragment(\@pos1,\@pos2);
sub Overlap_fragment{
	my $pos1_p = shift; 
	my $pos2_p = shift; 
	my $pos_p;  ## mix
	my $new_p; ## conjoin
	my $overlap_p; ## overlap
	
	my ($mix_size,$conjoin_size,$overlap_size);

	@$pos_p = (@$pos1_p,@$pos2_p);

	#make position in start-end order
	foreach my $p (@$pos_p) {
			($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
	}
	
	@$pos_p = sort {$a->[0] <=> $b->[0]} @$pos_p;
	my $mini_p = shift @$pos_p;
	push @$new_p, [ $mini_p->[0],$mini_p->[1] ];
	
	foreach my $p (@$pos_p) {
			if ( ($p->[0] - $new_p->[-1][1]) <= 0 ) { # conjoin
					if ($new_p->[-1][1] < $p->[1]) {
						push @$overlap_p, [$p->[0],$new_p->[-1][1]];
						$new_p->[-1][1] = $p->[1]; 
					}else{
						push @$overlap_p, [ $p->[0],$p->[1] ];
					}
					
			}else{  ## not conjoin
					push @$new_p, [ $p->[0],$p->[1] ];
			}
	}
	unshift @$pos_p,$mini_p;


	foreach my $p (@$pos_p) {
		$mix_size += $p->[1] - $p->[0] + 1;
	}
	foreach my $p (@$new_p) {
		$conjoin_size += $p->[1] - $p->[0] + 1;
	}
	foreach my $p (@$overlap_p) {
		$overlap_size += $p->[1] - $p->[0] + 1;
	}
	

	print STDERR "$mix_size\t$conjoin_size\t$overlap_size\n";
	return $overlap_p;
}


##check a seqeunce to be correct CDS-ORF or not
sub check_CDS{
	my $seq=shift;
	my $check=shift;
	$seq=~s/\s//g;
	$seq=uc($seq);
	my ($start,$end,$mid,$triple,$mask);
	$mid=1;
	my $len=length($seq);
	my $Ns=$seq=~tr/N//;
	$mask=$Ns/$len;
	$triple=1 if($len%3 == 0);
	$start=1 if($seq=~/^ATG/);
	$end=1 if($seq=~/TAA$|TAG$|TGA$/);
	for (my $i=3; $i<$len-3; $i+=3) {
		my $codon=substr($seq,$i,3);
		$mid=0 if($codon eq 'TGA' || $codon eq 'TAG' || $codon eq 'TAA');
	}
	#$start && $mid && $end && $triple 
	#!$start || !$end || !$mid || !$triple
	# !$start || !$end
	# !$mid
	# !$triple
	if ($check eq "wrong" && (!$start || !$end || !$mid || !$triple) ) {
		return 1;
	}elsif($check eq "right" && ($start && $mid && $end && $triple) ){
		return 1;
	}elsif($check eq "startend" && (!$start || !$end) ){
		return 1;
	}elsif($check eq "midstop" && !$mid ){
		return 1;
	}elsif($check eq "frameshift" && !$triple ){
		return 1;
	}else{
		return 0;
	}

}




1;

__END__