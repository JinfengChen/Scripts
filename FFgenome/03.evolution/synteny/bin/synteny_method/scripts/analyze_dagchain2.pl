#!/usr/bin/env perl

#analyze_dagchain.pl file.matchList file.aligncoords > analyze_nondag.out


use strict;

my ($match_file, $align_file) = @ARGV;

my ($matchlist, $ort1_2, $ord1, $ord2,) = get_orthologs($match_file);

my ($aln1_2, $daggene1, $dagsorted, $block) = get_dag_info($align_file);

print_header();

my (@data, %seen);

for my $line (@$matchlist){
    
    my ($ort1_ref, $ort2_ref) = ($line->[0], $line->[3]); 
    my ($ort1_gid, $ort2_gid) = ($line->[1], $line->[4]);
    my ($ort1_ord, $ort2_ord) = ($line->[2], $line->[5]);
    
    my $ort1_coord = $ort1_ref . ':' . $ort1_ord;
    my $ort2_coord = $ort2_ref . ':' . $ort2_ord;
    
    my $axis = $ort1_ref . '_' . $ort2_ref;                 
    my $block_id;
    
    if ($block->{$ort1_gid}->{$ort2_gid}){
	$block_id = join(';', @{$block->{$ort1_gid}->{$ort2_gid}} );
    } else {$block_id = 'none';}
    
    my ($up_aln1_gid, $up_aln1_ord, $dn_aln1_gid, $dn_aln1_ord)
         = get_up_down_dag($ort1_ref, $ort1_ord,);
    
    my $up_aln1_coord = $ort1_ref . ':' . $up_aln1_ord;
    my $dn_aln1_coord = $ort1_ref . ':' . $dn_aln1_ord;
    
    my ($up_dist1, $dn_dist1, $aln1_gap)
        = find_position($ort1_ref, $ort1_ord, 
	                $ort1_ref, $up_aln1_ord, 
			$ort1_ref, $dn_aln1_ord);
    
    my @up_aln2_gid;# = ('null');
    my @dn_aln2_gid;# = ('null');
    @up_aln2_gid = $aln1_2->{$up_aln1_gid} ? 
                   keys %{$aln1_2->{$up_aln1_gid}} :
		   ('null');

    @dn_aln2_gid = $aln1_2->{$dn_aln1_gid} ? 
                   keys %{$aln1_2->{$dn_aln1_gid}} :
		   ('null');

    for my $up_aln2_gid (@up_aln2_gid){
        my ($up_aln2_ref, $up_aln2_ord, $up_aln2_coord, $up_aln_block);

	unless ($up_aln2_gid eq 'null'){
            ($up_aln2_ref, $up_aln2_ord) = @{$ord2->{$up_aln2_gid}};

	    $up_aln_block 
	        = join(';', @{$block->{$up_aln1_gid}->{$up_aln2_gid}});

	}

        for my $dn_aln2_gid (@dn_aln2_gid){
	    my ($dn_aln2_ref, $dn_aln2_ord, $dn_aln2_coord, $dn_aln_block);
	    unless ($dn_aln2_gid eq 'null'){
	        ($dn_aln2_ref, $dn_aln2_ord) = @{$ord2->{$dn_aln2_gid}};

	        $dn_aln_block 
		    = join(';', @{$block->{$dn_aln1_gid}->{$dn_aln2_gid}});	    
	    }
	    
        my ($up_dist2, $dn_dist2, $aln2_gap)
            = find_position($ort2_ref, $ort2_ord,
		            $up_aln2_ref, $up_aln2_ord, 
		            $dn_aln2_ref, $dn_aln2_ord);
	
	my $aln = $block_id eq 'none' ? 0 : 1;

	my $data =
	join("\t",
             $ort1_gid,
             $ort1_ref,
             $ort1_ord,
             $ort2_gid,
             $ort2_ref,
             $ort2_ord,
	     $aln,
             $block_id,
             $up_dist2,
             $dn_dist2,
             $aln2_gap,
             $axis,     
             ), "\n";
	 
	 push @data, $data unless $seen{$data};  #remove redundancies
	 $seen{$data} = 1;
	 
	}
    } 
}

for my $line (@data){
    print $line, "\n";
}

sub find_position{
    my ($ort_ref, $ort_ord,
        $up_ref, $up_ord, 
        $dn_ref, $dn_ord) = @_;
    
    my ($gap, $up_dist, $dn_dist);
   
    $gap     = $up_ref  ne $dn_ref ? 'off' : abs($dn_ord - $up_ord);
    $up_dist = $ort_ref ne $up_ref ? 'off' : abs($ort_ord - $up_ord);
    $dn_dist = $ort_ref ne $dn_ref ? 'off' : abs($ort_ord - $dn_ord);

    return ($up_dist, $dn_dist, $gap);

}


sub get_up_down_dag {
    my ($ref, $ord, ) = @_; #$dagsorted) = @_;

    my @up = grep {$_->[1] < $ord } @{$dagsorted->{$ref}};
    my $up = $up[-1];
    my ($up_gid, $up_ord) = $up ? @$up : ('null', 'null');
    my @dn = grep {$_->[1] > $ord } @{$dagsorted->{$ref}};
    my $dn = $dn[0];

    my ($dn_gid, $dn_ord) = $dn ? @$dn : ('null', 'null');

    return ($up_gid, $up_ord, $dn_gid, $dn_ord);
}


sub get_orthologs {
    my $match_file = shift;
    my @lines;
    my (%ort1_2, 
	%ord1,
	%ord2);
    
    open my $IN, "<$match_file" or die "can't open $match_file\n";
    while(<$IN>){
        chomp;
	my ($ref1, $gid1, $ord1, undef, $ref2, $gid2, $ord2) = split /\t/;
        
	push @lines, [$ref1, $gid1, $ord1, $ref2, $gid2, $ord2];
	$ort1_2{$gid1}->{$gid2} = 1;
        $ord1{$gid1} = [$ref1, $ord1];
        $ord2{$gid2} = [$ref2, $ord2];
	
    }
    my @matchlist = sort {$a->[0] cmp $b->[0] || $a->[2] <=> $b->[2]} @lines;
    
    return (\@matchlist, \%ort1_2, \%ord1, \%ord2,);

}


sub get_dag_info {
    my $align_file = shift;
    my (%block, %aln1_2, %dagorder1, %daggene1,);
    my $block_id;
    
    open my $IN, "<$align_file" or die "can't open $align_file\n";
    while(<$IN>){
        chomp;
	split;
	if ($_[0] eq '##'){  #block header
	    my ($aln_num, $score, $pair_cnt, $refA, $refB) 
	        = ($_[-8], $_[-5], $_[-1], $_[2], $_[4]);
	    $aln_num =~ s/\#//;
	    $pair_cnt =~ s/\)://;
            my $strand = $_[5] eq '(reverse)' ? 'R' : 'F';
	    $block_id = $refB . '_' . $refA . '_' . 
	                $strand . '.' . $aln_num . 
			'(' . $pair_cnt . ')';  #eg: Zm1_Sb1_F.13(7)
	}
	else {
	    my ($ref1, $ref2) = ($_[0], $_[4]);
	    my ($gid1, $gid2) = ($_[1], $_[5]);
	    my ($ord1, $ord2) = ($_[2], $_[6]);
	    	    
	    push @{$block{$gid1}->{$gid2}}, $block_id;
	    
	    $aln1_2{$gid1}->{$gid2} = 1;
	    
	    push @{$dagorder1{$ref1}}, [$gid1, $ord1,] unless $daggene1{$gid1};
	    
	    $daggene1{$gid1} = [$ref1, $ord1,];
	}
    }
    
    my %dagsorted1 = sort_by_order(\%dagorder1);
    
    return (\%aln1_2, \%daggene1, \%dagsorted1, \%block);

}

sub sort_by_order {
    my ($dagorder) = @_;
    my %dagsorted;
    
    for my $ref (keys %$dagorder){
	my @sorted = sort { $a->[1]<=>$b->[1] } @{$dagorder->{$ref}};
	push @{$dagsorted{$ref}}, @sorted;
    }
    return %dagsorted;
}


sub print_header {
    print
    join
        ("\t",
        'ort1_gid',
	'ort1_ref',
	'ort1_ord',
        'ort2_gid',
	'ort2_ref',
	'ort2_ord',
	'align?',
        'block',
	'up2_dist',
	'dn2_dist',
	'gap2',
        'axis',
        ),"\n";
}
