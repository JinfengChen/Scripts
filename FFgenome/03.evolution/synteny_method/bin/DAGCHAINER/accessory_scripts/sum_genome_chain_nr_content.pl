#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Overlap_piler;


my $usage = "usage: $0 chain_spans_file [verbose]\n\n";

my $chain_coords_file = $ARGV[0] or die $usage;
my $verbose_flag = $ARGV[1] || 0;

main: {

	my %genomeA_data;
	my %genomeB_data;

	open (my $fh, $chain_coords_file) or die "Error, cannot open file $chain_coords_file";
	while (<$fh>) {
		if (/^\#/ || ! /\w/) {
			die "Error, cannot parse chain span line: $_";
		}
		chomp;
		my ($molA, $lendA, $rendA, $molB, $lendB, $rendB, $orient, $num_genes) = split (/\t/);
		
		push (@{$genomeA_data{$molA}}, [$lendA, $rendA]);
		push (@{$genomeB_data{$molB}}, [$lendB, $rendB]);

	}
	close $fh;

	my $sum_bp_A = &sum_coord_piles("A", \%genomeA_data);
	my $sum_bp_B = &sum_coord_piles("B", \%genomeB_data);

	print "GenomeA: $sum_bp_A\n"
		. "GenomeB: $sum_bp_B\n";

	exit(0);

}


####
sub sum_coord_piles {
	my ($org_token, $data_href) = @_;

	my $sum_bp = 0;
	
	foreach my $mol (keys %$data_href) {
		
		my @coordpairs = @{$data_href->{$mol}};
		my @piles = &Overlap_piler::simple_coordsets_collapser(@coordpairs);

		foreach my $pile (@piles) {
			my ($lend, $rend) = @$pile;
			my $pile_len = $rend - $lend + 1;
			$sum_bp += $pile_len;
			if ($verbose_flag) {
				print "$org_token\t$mol\t$lend\t$rend\n";
			}
		}
	}
	
	return ($sum_bp);
}
	
		
