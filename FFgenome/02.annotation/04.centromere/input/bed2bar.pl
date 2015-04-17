#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for window view read hit counts on genome.
Alignment of bed format should be sorted using msort (msort -k 1,n2 DNAmethly.bed > DNAmethly.bed.sorted).
=cut

use Getopt::Long;

my %opt;
GetOptions(\%opt,"window:s","bar:s","bed:s","sorting","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -win 50 -bar DNAmethly.bar.txt -bed DNAmethly.bed.sorted\n";
   exit();
}
$opt{window} ||= 50;

if ($opt{sorting}){
   system ("msort -k 1,n2 $opt{bed} > $opt{bed}.sorted");
   system ("rm $opt{bed}");
   system ("mv $opt{bed}.sorted $opt{bed}");
}
my %hash;
my $chr="chr1";
open IN, "$opt{bed}" or die "$!";
open OUT, ">$opt{bar}" or die "$!";

