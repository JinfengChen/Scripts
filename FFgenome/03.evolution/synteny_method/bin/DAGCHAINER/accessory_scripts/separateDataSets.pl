#!/usr/local/bin/perl

use strict;

my $filename = $ARGV[0] or die "usage: $0 matchFile";

my %files;

open LADEDA, $filename or die $!;
while (<LADEDA>) {
    my @x = split (/\t/);
    if (! (scalar(@x) >= 8)) {  next;}
    my ($mol1, $mol2) = ($x[0], $x[4]);
    
    my @mols = sort ($mol1, $mol2);
    my $molPairKey = join ("vs", @mols);
    
    my $filename = "$filename.$molPairKey";
    unless ($files{$filename}) {
	if (-s $filename) {
	    unlink $filename;
	}
	$files{$filename} = 1;
    }
    open (FILE, ">>$filename") or die $!;
    print FILE $_;
    close FILE;
}

close LADEDA;
