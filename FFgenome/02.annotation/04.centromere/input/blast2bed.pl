#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blastm8:s","help");


my $help=<<USAGE;
perl $0 --blastm8 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 


my $head=$1 if ($opt{blastm8}=~/(.*)\.blast/);
open OUT, ">$head.bed" or die "$!";
open IN, "$opt{blastm8}" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    my $strand= $unit[9]-$unit[8] > 0 ? "+" : "-";
    if ($unit[9]-$unit[8] > 0){
       print OUT "$unit[1]\t$unit[8]\t$unit[9]\t$unit[2]\t$unit[10]\t$strand\n";
    }else{
       print OUT "$unit[1]\t$unit[9]\t$unit[8]\t$unit[2]\t$unit[10]\t$strand\n";
    }
}
close IN;
close OUT;
