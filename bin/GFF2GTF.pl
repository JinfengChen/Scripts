#!/usr/bin/perl

use warnings;
use strict;
while (<>) {
chomp;
my @parts = split /\t/;
if ( !(($parts[2] eq 'three_prime_UTR')|($parts[2] eq 'five_prime_UTR')|($parts[2] eq 'exon')|($parts[2] eq 'CDS')) ) {
next;
}
elsif ( $parts[2] eq 'three_prime_UTR') {
$parts[2] = '3UTR'; 
}
elsif ($parts[2] eq 'five_prime_UTR') {
$parts[2] = '5UTR';
}
$parts[8] =~ s/.*Parent=((\w+)\.\w*).*/gene_id \"$2\"; transcript_id \"$1\";/sg;
$_ = join "\t" , @parts;
print "$_\n";
}


