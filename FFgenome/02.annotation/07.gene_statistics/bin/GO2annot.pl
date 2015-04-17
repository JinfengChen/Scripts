#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Convert BGI *.pep.iprscan.gene.GO to blast2go *.annot.
The latter can be used to do GO enrichment analysis by online blast2go.
Run: perl GO2annot.pl *.pep.iprscan.gene.GO
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 

open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].annot" or die "$!";
while (<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    my $head=shift @unit;
    my $num=shift @unit;
    foreach (@unit){
        my @go=split(";",$_);
        $go[1]=~s/^\s+//g;
        print OUT "$head\t$go[0]\t$go[1]\n";
    }
}
close IN;
close OUT;
