#!/usr/local/bin/perl

use strict;

my $usage = "usage: matchFile diagFile\n";

my $matchFile = $ARGV[0] or die $usage;
my $diagFile = $ARGV[1] or die $usage;

my %chromoPairs;
open (DIAG, "$diagFile") or die "Cannot open $diagFile\n";
while (<DIAG>) {
    my @x = split (/\t/);
    my ($chrA, $chrB) = sort ($x[0], $x[4]);
    $chromoPairs{"$chrA,$chrB"} = 1;
}
close DIAG;

open (MATCHES, "$matchFile") or die "Cannot open $matchFile\n";
while (<MATCHES>) {
    my @x = split (/\t/);
    my ($chrA, $chrB) = sort ($x[0], $x[4]);
    if ($chromoPairs{"$chrA,$chrB"}) {
	print;
    }
}
close MATCHES;
