#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","feature:s","bed:s","help");


my $help=<<USAGE;
Convert GFF3 format into BED format with user specified feature.
-gff:  GFF3 format file, input
-bed:  BED format file, output
-feature: Feature to be dealed with

Run: perl GFF2BED.pl -gff FF.mRNA.gff -feature mRNA -bed FF.mRNA.bed

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

open IN, "$opt{gff}" or die "$!";
open OUT, ">$opt{bed}.unsort" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[2] eq $opt{feature} and ($unit[8] =~/ID=(.*?);/ or $unit[8] =~/Parent=(.*?);/)){
        my $start=$unit[3]-1; ### GFF start is 1-based, but BED start is 0-based. See BEDtools Manual for detail.
        print OUT "$unit[0]\t$start\t$unit[4]\t$1\t$unit[5]\t$unit[6]\n"; 
    }
}
close IN;
close OUT;

system ("msort -k 1,n2 $opt{bed}.unsort > $opt{bed}");
system ("rm $opt{bed}.unsort");
 
