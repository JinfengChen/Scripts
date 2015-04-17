#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Std;
use Storable qw(dclone);

my $DEBUG = 0;
$|=1;

use vars qw ($opt_K $opt_i $opt_H $opt_G $opt_T $opt_M $opt_F $opt_D $opt_S $opt_m $opt_A $opt_E $opt_v $opt_h $opt_c $opt_d $opt_s $opt_C $opt_S $opt_Z $opt_x $opt_g $opt_e $opt_o $opt_I $opt_g);

&getopts ('Ic:d:HG:TM:F:D:S:m:A:E:vhsi:CS:Z:x:o:e:g:K');

my $usage =  <<_EOH_;

############################# Options ###############################
#

## Required:
# -i input file
#          input file has format:
#          molecule_1 <tab> accession_1 <tab> end5_1 <tab> end3_1 <tab> molecule_2 <tab> accession_2 <tab> end5_2 <tab> end3_2 <tab> P-value


## DAG chaining scoring parameters:
# -o gap open penalty (default: -0f)
# -e gap extension penalty (default: -3f)
# -g length of a gap in bp (default: 10000)  (avg distance expected between two syntenic genes) 
# -M Maximum match score (default: 50) otherwise, -log(evalue)
#     -Z define constant match score ** use in place of -M
# -D maximum distance allowed between two matches in basepairs. (default: 200000)

## Input data filtering:
# -E Maximum E-value (default 1e-5)


## Output filtering:
# -x Minimum alignment score (alignment pursued until scores fall below this range.)  
#    default:  MIN_ALIGNMENT_SCORE = (int) (MIN_ALIGN_LEN * 2.5 * -GAP_PENALTY)
# -I ignore tandem duplication alignments (overlapping, same mol alignments) (requires -s).
# -A Minium number of Aligned Pairs (default: 6)
# -T only Tandem alignments (implies -s)


## Include/Exclude certain molecule comparisons:
# -s include self comparisons.


## Others:
# -v verbose
# -h print this option menu and quit
###################### Process Args and Options #####################

_EOH_

    ;

die $usage if $opt_h;

my $inputFile = $opt_i or die $usage;
my $SEE = $opt_v;
my $GAP_LENGTH = $opt_g || 10000;
my $GAP_OPEN_PENALTY = (defined($opt_o)) ? $opt_o : 0;
my $GAP_EXTENSION_PENALTY = (defined($opt_e)) ? $opt_e : -3;
my $MAX_MATCH_SCORE = $opt_M || 50;
my $MAX_DIST_BETWEEN_MATCHES = (defined ($opt_D)) ? $opt_D : 200000;

my $INCLUDE_SELF = $opt_s;
my $IGNORE_TANDEM = $opt_I;
my $TANDEM_ONLY = $opt_T;
if ($TANDEM_ONLY) {
    $INCLUDE_SELF = 1;
}
if ($IGNORE_TANDEM && $TANDEM_ONLY) {
    die "Can't look for tandems and ignore them at the same time.n";
}
my $CONSTANT_MATCH_SCORE = (defined $opt_Z) ? $opt_Z : undef;
if ($CONSTANT_MATCH_SCORE) {
    $MAX_MATCH_SCORE = $CONSTANT_MATCH_SCORE; # for auto-determine min alignment score cutoff
}

our $KEEP_DELCHER_FILES = $opt_K; #hidden opt for troubleshooting.

my $uname = $ENV{HOSTTYPE};

my $prog = $0;
my $progpath;
if ($prog =~ /\//) {
    $progpath = $prog;
    $progpath =~ s/[^\/]+$//;
}
print "PROGPATH: $progpath\n";

my $MIN_NUM_ALIGNED_PAIRS = $opt_A || 6;
my $max_e_value = defined($opt_E) ? $opt_E : 1e-05;
my $MIN_ALIGNMENT_SCORE = $opt_x;
unless ($MIN_ALIGNMENT_SCORE) {
    $MIN_ALIGNMENT_SCORE = int ($MIN_NUM_ALIGNED_PAIRS * 0.5 * $MAX_MATCH_SCORE);
}

## Primary Global vars in use:
my %mol_pair_to_list; ## molPairKey -> aref of matches
my %acc_info; ## accession -> featureStruct
my %acc_pair_to_match; ## accPairKey -> matchStruct



main: {
	 
	## Parse input match data, populate data structures:
	$|++;
	print STDOUT "Parsing inputFile\n";
	&parse_input_file($inputFile);
	
    ## open output file and force buffer flushing.
    open (ALIGNCOORDS, ">$inputFile.aligncoords") or die "Can't open $inputFile.aligncoords";
    my $ref = select ALIGNCOORDS;
    $|++;
    select $ref;
    print STDOUT "Done parsing inputfile\n";
	
    
    ## Perform the DAG chaining for each molecule pair (do +/- orientations separately)
    foreach my $molpair (sort keys %mol_pair_to_list) {
		print "molpair: $molpair\n" if $DEBUG||$SEE;
		my ($mol_1, $mol_2) = split (/,/, $molpair);
		print "***** Comparing $mol_1 to $mol_2 *****\n";
		my $match_list_aref = $mol_pair_to_list{$molpair};
		
		## Create input file to cpp program.
		my $filename = "$mol_1.vs.$mol_2.delcher..input";
		
		my @pairIndexToAccs;
		open (ART, ">$filename") or die "Can't open file.\n";
		my $pairID = 0;
		foreach my $match (@$match_list_aref) {
			my ($feature_A, $feature_B, $e_value) = ( $match->{feature_A},
													  $match->{feature_B},
													  $match->{e_value} 
													  );
			my $midpt_A = $feature_A->{mid};
			my $midpt_B = $feature_B->{mid};
			
			my $e_value = $match->{e_value};
			
			my $score = scoringF($e_value);
			
			printf ART ("%d\t%d\t%d\t%f\n", $pairID, $midpt_A, $midpt_B, $score);
			
			$pairIndexToAccs[$pairID] = [$feature_A->{acc}, $feature_B->{acc}];
			$pairID++;
		}
		
		close ART;
		
		&run_DAG_chainer($mol_1, $mol_2, $filename, \@pairIndexToAccs, ""); ## forward direction
		
		
		&run_DAG_chainer($mol_1, $mol_2, $filename, \@pairIndexToAccs, "-r"); ## revcomp mol2
		
		unlink ($filename) unless ($KEEP_DELCHER_FILES); #remove tempfile.
		
    } #end of molecule.
    
    close ALIGNCOORDS;
}


####
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

####
sub print_alignment {
    my ($mol1, $mol2, $match_header, $align_list_aref) = @_;
    my $ignore_alignment = 0;
    my $ignore_reason = "";
	
    my $num_aligned_pairs = scalar (@$align_list_aref);
    print "# $mol1 vs. $mol2 $match_header $num_aligned_pairs aligned pairs.\n";
    if ($num_aligned_pairs < $MIN_NUM_ALIGNED_PAIRS) { #more than one pair aligned.
		$ignore_alignment = 1;
		$ignore_reason = " ignoring alignment.  Only $num_aligned_pairs, threshold=$MIN_NUM_ALIGNED_PAIRS";
    }
    
    my $alignment_text = "## alignment $mol1 vs. $mol2 $match_header (num aligned pairs: $num_aligned_pairs):\n";
    my $IS_TANDEM = 0;
    my $IS_OVERLAPPING = 0;
    if ((!$ignore_alignment) && ($mol1 eq $mol2) && ($IGNORE_TANDEM || $TANDEM_ONLY)) {
		my @acc1_coords;
		my @acc2_coords;
		print "## checking to see if should ignore a tandem aligment.\n" if $SEE;
		## walk thru and see if any accession is showing up multiple times in alignment:
		my %seen;
		foreach my $alignedPair (@$align_list_aref) {
			my ($mol1_acc, $mol2_acc) = ($alignedPair->{acc_A}, $alignedPair->{acc_B});
			print "$mol1_acc,$mol2_acc\n" if $SEE;
			if (($mol1_acc ne '-' && $seen{$mol1_acc}) || ($mol2_acc ne '-' && $seen{$mol2_acc})) {
				print "SEEN!!, ignoring.\n" if $SEE;
				$IS_TANDEM = 1;
				last;
			} else {
				$seen{$mol1_acc} = 1; 
				$seen{$mol2_acc} = 1;
			}
			## get coords for is_overlap determination.
			my $accPairKey =  &get_acc_pair_key($mol1_acc, $mol2_acc);
			my $match = $acc_pair_to_match{$accPairKey};
			my ($feature_A, $feature_B) = ($match->{feature_A}, $match->{feature_B});
			if ($mol1_acc eq $feature_B->{acc}) {
				($feature_A, $feature_B) = ($feature_B, $feature_A); #swap 'em
			}
			
			my ($a_end5, $a_end3, $b_end5, $b_end3) = ($feature_A->{end5}, $feature_A->{end3}, $feature_B->{end5}, $feature_B->{end3});
			push (@acc1_coords, $a_end5, $a_end3);
			push (@acc2_coords, $b_end5, $b_end3);
			
		}
		
		@acc1_coords = sort {$a<=>$b} @acc1_coords;
		my $lend_A = shift @acc1_coords;
		my $rend_A = pop @acc1_coords;
		
		@acc2_coords = sort {$a<=>$b} @acc2_coords;
		my $lend_B = shift @acc2_coords;
		my $rend_B = pop @acc2_coords;
		
		if ($lend_A <= $rend_B && $rend_A >= $lend_B) { #overlap
			$IS_OVERLAPPING = 1;
		}
		
    }
    
    if ($IGNORE_TANDEM && ($IS_TANDEM || $IS_OVERLAPPING)) {
		$ignore_alignment = 1;
		$ignore_reason = "ignoring tandem or overlapping self alignments.";
    } elsif ($TANDEM_ONLY && !$IS_TANDEM) {
		$ignore_alignment = 1;
		$ignore_reason = "ignoring, it's not tandem and only want tandem.";
    }
    
    if ($ignore_alignment) {
		print "# $ignore_reason\n";
    } else {
		print ALIGNCOORDS $alignment_text;
		print $alignment_text if $SEE;
		
		foreach my $alignedPair (@$align_list_aref) {
			my ($acc_A, $acc_B, $dag_position_score) = ($alignedPair->{acc_A},
														$alignedPair->{acc_B},
														$alignedPair->{dag_position_score}
														);
			
			my $accPairKey =  &get_acc_pair_key($acc_A, $acc_B);
			my $match = $acc_pair_to_match{$accPairKey};
			
			my ($feature_A, $feature_B, $e_value) = ($match->{feature_A},
													 $match->{feature_B},
													 $match->{e_value}
													 );
			
			
			my ($mol_1, $acc_1, $end5_1, $end3_1) = ($feature_A->{mol},
													 $feature_A->{acc},
													 $feature_A->{end5},
													 $feature_A->{end3}
													 );
			
			
			
			
			
			my ($mol_2, $acc_2, $end5_2, $end3_2) = ($feature_B->{mol},
													 $feature_B->{acc},
													 $feature_B->{end5},
													 $feature_B->{end3});
			
			my $outline = sprintf ("%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%e\t%d\n", 
								   $mol1, $acc_1, 
								   $end5_1, $end3_1, 
								   $mol2, $acc_2, 
								   $end5_2, $end3_2, 
								   $e_value, $dag_position_score);
			
			print ALIGNCOORDS $outline;
			print $outline if $SEE;
		}
    }
}


####
sub get_mol_pair_aref {
    my @mols  = @_;
    @mols = sort @mols;
    my $key = join (",", @mols);
    if (my $aref = $mol_pair_to_list{$key}) {
		return ($aref);
    } else {
		print "Getting molpair: @mols\n" if $SEE;
		my $aref = [];
		$mol_pair_to_list{$key} = $aref;
		return ($aref);
    }
}

####
sub count_num_pairs {
    my $aref = shift;
    my $count = 0;
    foreach my $pair (@$aref) {
		my ($a, $b) = @$pair;
		#print "a: $a\tb: $b\n";
		unless ($a eq "-" || $b eq "-") {
			$count++;
		}
    }
    return ($count);
}


####
sub get_acc_pair_key {
    my (@accs) = @_;
    @accs = sort (@accs);
    my $key = join (",", @accs);
    return ($key);
}


####
sub parse_input_file {
    my ($file) = shift;
    
    my %seenAcc; #track accessions examined already.
    open (FILE, $file) or die "Can't open $file";
    while (<FILE>) {
		chomp;
		unless (/\w/) { next;}
		my @x = split (/\t/);
		my ($mol_1, $acc_1, $end5_1, $end3_1, $mol_2, $acc_2, $end5_2, $end3_2, $e_value) = @x;
		
		#unless ( ($mol_1 eq "10_319" && $mol_2 eq "18_68") 
		#	 || 
		#	 ($mol_2 eq "10_319" && $mol_1 eq "18_68") 
		#	 ) { next;}
		
		if ($e_value < 1e-250) {
			## can't take logs of 0
			$e_value = 1e-250;
		}
		
		## run a few checks, see if we want this record.
		if ($acc_1 eq $acc_2) { next;} #no self comparisons.
		unless ($e_value <= $max_e_value) { next;}
		if ($mol_1 eq $mol_2 && !$INCLUDE_SELF) {next;}
		if ($mol_1 ne $mol_2 && $TANDEM_ONLY) { next;}
		
		my ($feature_1, $feature_2);
		if (! $seenAcc{$acc_1}) {
			$feature_1 = &store_acc_info($mol_1, $acc_1, $end5_1, $end3_1);
			$seenAcc{$acc_1} = 1;
		} else {
			$feature_1 = $acc_info{$acc_1};
		}
		
		if (! $seenAcc{$acc_2}) {
			$feature_2 = &store_acc_info($mol_2, $acc_2, $end5_2, $end3_2) unless $seenAcc{$acc_2};
			$seenAcc{$acc_2} = 1;
		} else {
			$feature_2 = $acc_info{$acc_2};
		}
		
		## keep match data in order by molecule name (lexically).
		my ($feature_A, $feature_B, $mol_A, $mol_B);
		if ($mol_1 lt $mol_2) {
			($feature_A, $feature_B) = ($feature_1, $feature_2);
			($mol_A, $mol_B) = ($mol_1, $mol_2);
		} else {
			($feature_A, $feature_B) = ($feature_2, $feature_1);
			($mol_A, $mol_B) = ($mol_2, $mol_1);
		}
		
		## If on the same molecule, order according to midpoint to keep all data points in the same side (mirror image) of the dot plot.
		if ($mol_A eq $mol_B && $feature_A->{mid} > $feature_B->{mid}) {
			($feature_A, $feature_B) = ($feature_B, $feature_A);
		}
		
		my $accPairKey =  &get_acc_pair_key($acc_1, $acc_2);
		my $match = $acc_pair_to_match{$accPairKey};
		if (ref $match) {
			if ($match->{e_value} > $e_value) {
				$match->{e_value} = $e_value; #take the lowest e_value for the match pair.  Important when there are multiple HSPs reported between two accessions.
			}
		} else {
			$match = {  feature_A => $feature_A,
						feature_B => $feature_B,
						e_value => $e_value 
						};
			
			## Store match based on molecule-pair comparison:
			my $mol_pair_aref = &get_mol_pair_aref($mol_A, $mol_B);
			push (@$mol_pair_aref, $match);
			$acc_pair_to_match{$accPairKey} = $match;
			
		}
    }
    close FILE;
	
	
}

####
sub run_DAG_chainer {
    my ($mol1, $mol2, $filename, $pairIndexToAccs_aref, $reverseOrientFlag) = @_;
	
    ## Run cpp program.
    
	if ($SEE) {
		print "\n\n\nDAGCHAINERcpp input:\n";
		system "cat $filename";
	}
	
    my $tmpFile = ".$$.tmpOut";
    my $cmd = "${progpath}dagchainer.$uname -G $GAP_LENGTH -O $GAP_OPEN_PENALTY -E $GAP_EXTENSION_PENALTY -S $MIN_ALIGNMENT_SCORE -D $MAX_DIST_BETWEEN_MATCHES  -F $filename $reverseOrientFlag > $tmpFile";
    print "CMD: $cmd \n";
    my $ret = system $cmd;
    if ($ret) {
		die "ERROR, couldn't run command\n\n$cmd\n(ret: $ret)\n\nDid you recompile ${progpath}dagchainer.$uname for your OS?\n\n";
    }
    my $all_alignments;
    open (OUT, "<$tmpFile");
    while (<OUT>) {
		print if $SEE;
		$all_alignments .= $_;
    }
    close OUT;
    unlink $tmpFile;
    
    
    my @individual_alignments = split (/>/, $all_alignments);
    shift @individual_alignments; #rid header that lacks alignment
    foreach my $individual_alignment (@individual_alignments) {
		my @align;
		my @matches = (split (/\n/, $individual_alignment));
		
		my $match_header = shift @matches;
		
		foreach my $match (@matches) {
			$match =~ s/^\s+//;
			
			if ($match =~ /^\d+:/) {
				chomp $match;
				my @x = split (/\s+/, $match);
				my ($index, $pairID, $pos1, $pos2, $match_score, $dag_chain_score) = @x;
				my $acc_pair_aref = $pairIndexToAccs_aref->[$pairID];
				my ($accA, $accB) = @$acc_pair_aref;
				push (@align, {acc_A=> $accA,
							   acc_B=> $accB,
							   dag_position_score => $dag_chain_score
							   }
					  );
			}
		}
		
		if ($reverseOrientFlag) {
			$match_header = "(reverse) $match_header";
		}
		
		&print_alignment($mol1, $mol2, $match_header, \@align);
    }
}



####
sub scoringF {
    my $evalue = shift;
    if ($CONSTANT_MATCH_SCORE) {
		return ($CONSTANT_MATCH_SCORE);
    } else {
		my $matchScore = -log10($evalue);
		$matchScore *=10;
		$matchScore = int($matchScore +.5);
		$matchScore /= 10;
		$matchScore = $MAX_MATCH_SCORE if $matchScore > $MAX_MATCH_SCORE;
		return ($matchScore);
    }
}


####
sub store_acc_info {
    my ($molecule, $accession, $end5, $end3) = @_;
    
    my $midPt = int (($end5+$end3)/2 + 0.5);
	
    my $struct = {mol=>$molecule,
				  acc=>$accession, 
				  end5=>$end5,
				  end3=>$end3,
				  mid=>$midPt
				  }; 
    
    ## map entry to accession:
    $acc_info{$accession} = $struct;
	
    return ($struct);
	
}
