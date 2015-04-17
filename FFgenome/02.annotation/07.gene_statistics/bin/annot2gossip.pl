#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Convert blast2go *.annot to gossip *.id.
The latter can be used to do GO enrichment analysis by gossip.

USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 

my %hash;
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].id" or die "$!";
print OUT "GoStat IDs Format Version 1.0\n";
while (<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    my $goid;
    if ($unit[1]=~/GO:0*(\d+)/){
       $goid=$1;
    }
    if (exists $hash{$unit[0]}){
        $hash{$unit[0]}.=",$goid";
    }else{
        $hash{$unit[0]}=$goid;
    }
}

foreach (sort keys %hash){
    print OUT "$_\t$hash{$_}\n";

}
close IN;
close OUT;
