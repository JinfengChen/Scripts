#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
	chomp;
	my @x = split (/\t/);
	my ($contigA, $lendA, $rendA, $contigB, $lendB, $rendB, @res) = @x;

	my $lenA = $rendA - $lendA + 1;
	my $lenB = $rendB - $lendB + 1;

	my $ratio = $lenA / $lenB;
	
	print "$contigA\t$lenA\t$contigB\t$lenB\t$ratio\n";
}

exit(0);

