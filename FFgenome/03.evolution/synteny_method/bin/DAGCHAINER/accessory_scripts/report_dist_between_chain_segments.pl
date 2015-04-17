#!/usr/bin/env perl

use strict;
use warnings;


my $prev_struct = undef;

while (<>) {
	chomp;
	if (/^\#/) {
		$prev_struct = undef;
		next;
	}

	my ($contig_A, $acc_A, $lend_A, $rend_A, $contig_B, $acc_B, $lend_B, $rend_B, $Evalue, $chain_score) = split (/\t/);
	
	my $curr_struct = {
		mid_A => ($lend_A + $rend_A) / 2,
		
		mid_B => ($lend_B + $rend_B) / 2,
	};
	
	if ($prev_struct) {
		## compare them, report deltas between them.
		
		my $delta_A = abs($curr_struct->{mid_A} - $prev_struct->{mid_A});
		
		my $delta_B = abs($curr_struct->{mid_B} - $prev_struct->{mid_B});

		print "$contig_A\t$delta_A\t$contig_B\t$delta_B\n";
	}
	
	$prev_struct = $curr_struct;
}

exit(0);

