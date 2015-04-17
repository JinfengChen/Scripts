#!/usr/local/bin/perl

use strict;

my $usage = "usage: $0 windowLinkSize < inputFile\n\n"
    . "In all pairwise comparisons between chromosomes, multiple matches to the same accession are linked together if they fall within the windowLinkSize specified.  Only the match with the smallest E-value within each linked match set will pass thru the filter.\n\n";

my $windowLinkSize = $ARGV[0] or die $usage;

my %chromoPairToMatchList;

while (<STDIN>) {
    my $inputLine = $_;
    chomp;
    my @x = split (/\t/);
    my ($chrA, $accA, $posA_end5, $posA_end3, $chrB, $accB, $posB_end5, $posB_end3, $evalue) = @x;
    
    my $featA = { chr => $chrA,
		  acc => $accA,
		  mid => int (($posA_end5 + $posA_end3)/2 + 0.5)
		  };
    my $featB = { chr => $chrB,
		  acc => $accB,
		  mid => int (($posB_end5 + $posB_end3)/2 + 0.5)
		  };
    
    my ($chrPair, $match);
    if ($chrA gt $chrB) {
	($featA, $featB) = ($featB, $featA);
	$chrPair = "$chrB,$chrA";
    } else {
	$chrPair = "$chrA,$chrB";
    }
    my $match = { featA => $featA,
		  featB => $featB,
		  input => $inputLine,
		  evalue => $evalue
		  };
    my $list_aref = $chromoPairToMatchList{$chrPair};
    unless (ref $list_aref eq "ARRAY") {
	$list_aref = $chromoPairToMatchList{$chrPair} = [];
    }

    push (@$list_aref, $match);
}



## Examine each chromosome pairwise comparison, link matches
foreach my $chromoPair (sort keys %chromoPairToMatchList) {
    my $matchList_aref = $chromoPairToMatchList{$chromoPair};
    
    my @filteredMatches = &filter_matches ("featA", "featB", $matchList_aref);

    @filteredMatches = &filter_matches ("featB", "featA", \@filteredMatches);
    
    
    foreach my $match (@filteredMatches) {
	print $match->{input};
    }
    
}


####
sub filter_matches {
    my ($accKey, $coordKey, $list_aref) = @_;
    
    @$list_aref = sort {
	$a->{$accKey}->{acc} cmp $b->{$accKey}->{acc}
	||
	    $a->{$coordKey}->{mid} <=> $b->{$coordKey}->{mid}
    } @$list_aref;
    
        
    my @matchBins = ( [ $list_aref->[0] ] );
    foreach (my $i=1; $i <= $#$list_aref; $i++) {
	my $currMatch = $list_aref->[$i];

	my $prevMatchList_aref = $matchBins[$#matchBins];
	my $prevMatch = $prevMatchList_aref->[$#$prevMatchList_aref];
	my $prevAcc = $prevMatch->{$accKey}->{acc};
	my $prevCoordB = $prevMatch->{$coordKey}->{mid};
	
	my $currAcc = $currMatch->{$accKey}->{acc};
	my $currCoordB = $currMatch->{$coordKey}->{mid};
	
	if ( ($currAcc ne $prevAcc) || 
	     ($currAcc eq $prevAcc && ($currCoordB - $prevCoordB > $windowLinkSize) )
	     ) {
	    push (@matchBins, [$currMatch]); # start a new bin
	} else {
	    push (@$prevMatchList_aref, $currMatch); # add to existing bin
	}
    }
    
    
    ## iterate thru match bins, only print the entry with the lowest evalue
    
    my @retMatches;
    foreach my $bin (@matchBins) {
	my $binSize = scalar (@$bin);
	if ($binSize == 1) {
	    # only one entry, print it:
	    push (@retMatches, $bin->[0]);
	} else {
	    ## sort by evalue, print lowest evalue entry
	    @$bin = sort {$a->{evalue} <=> $b->{evalue}} @$bin;
	    push (@retMatches, $bin->[0]);
	}
    }

    return (@retMatches);
}


exit(0);


		     
		  
