#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","project:s","help");


my $help=<<USAGE;

Convert bed format file to mcscan bed format.
Example mcscan gff:
Os01 Os01g0100100 1903 9816
Sb01 Sb01g001001  2000 3000

Run: perl $0 -bed Rice.bed -project Os
-bed:
-project: prefix for chromosome name in mcscan gff, Os01

USAGE

print $help and exit if (keys %opt < 1); 

open IN, "$opt{bed}" or die "$!";
open OUT, ">$opt{project}.bed" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    $unit[0]=~s/[a-zA-Z]//g;
    $unit[0]=~s/\_//g;
    $unit[0] = "0".$unit[0] if (length $unit[0] == 1);
    $unit[0] =$opt{project}.$unit[0];
    #print OUT "$unit[0]\t$unit[3]\t$unit[1]\t$unit[2]\n";
    print OUT "$unit[0]\t$unit[1]\t$unit[2]\t$unit[3]\n";
    #my $line=join("\t",@unit);
    #print OUT "$line\n";
}
close IN;
close OUT;
