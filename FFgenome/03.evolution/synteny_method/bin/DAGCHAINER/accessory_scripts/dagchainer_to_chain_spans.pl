#!/usr/bin/env perl

use strict;
use warnings;


my @chains;
my $curr_chain = undef;

while (<>) {
	chomp;
	my $line = $_;
	if (/^\#/) {
		
		my @x = split (/\s+/, $line);
		my $accA = $x[2];
		my $accB = $x[4];
		my $orientation = ($line =~ /\(reverse\)/) ? '-' : '+';
		$line =~ /num aligned pairs: (\d+)/ or die "Error, cannot find num aligned pairs";
		my $num_aligned_pairs = $1;
		
		$curr_chain = {
			accA => $accA,
			accB => $accB,
			orient => $orientation,
			num_gene_pairs => $num_aligned_pairs,
			coordsA => [],
			coordsB => [],
			
		};
		
		push (@chains, $curr_chain);
		
	}
	else {
		my ($contigA, $geneA, $lendA, $rendA, 
			$contigB, $geneB, $lendB, $rendB,
			$evalue, $dagchainscore) = split (/\t/);

		push (@{$curr_chain->{coordsA}}, $lendA, $rendA);
		push (@{$curr_chain->{coordsB}}, $lendB, $rendB);

	}
	
}

foreach my $chain (@chains) {
	
	my $accA = $chain->{accA};
	my $accB = $chain->{accB};
	my $orient = $chain->{orient};
	my $num_gene_pairs = $chain->{num_gene_pairs};

	my @coordsA = @{$chain->{coordsA}};
	@coordsA = sort {$a<=>$b} @coordsA;
	my $lendA = shift @coordsA;
	my $rendA = pop @coordsA;
	

	my @coordsB = @{$chain->{coordsB}};
	@coordsB = sort {$a<=>$b} @coordsB;
	my $lendB = shift @coordsB;
	my $rendB = pop @coordsB;

	print "$accA\t$lendA\t$rendA\t$accB\t$lendB\t$rendB\t$orient\t$num_gene_pairs\n";

}


exit(0);

	
	
	
