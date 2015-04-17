#!/usr/bin/perl

my $help=<<USAGE;
This script is designed to create config file for drawRegion.pl.
Read CompareRegion.txt, output species, fregment, length to config file.

Example CompareRegion.txt:
OB_Scaffold000002_1830930_2143221       Os_chr09_18409023_18805233

Example drawRegion.config:
Oryza brachyantha       OB_Scaffold000002_1830930_2143221       312292
Oryza sativa    Os_chr09_18409023_18805233      396211

Run: perl config.pl CompareRegion.txt

USAGE

if (@ARGV < 1){
   print "$help\n";
   exit();
}

my %species;
$species{OB}="Oryza brachyantha";
$species{Os}="Oryza sativa";

open OUT, ">drawRegion.config" or die "$!";
open IN, "$ARGV[0]" or die "$!";
    while(<IN>){
       chomp $_;
       my @unit=split("\t",$_);
       foreach (@unit){
          my @temp=split("_",$_);
          my $len=$temp[3]-$temp[2]+1;
          my $spe=$species{$temp[0]};
          print OUT "$spe\t$_\t$len\n";
       }        
    }
close IN;
close OUT;
