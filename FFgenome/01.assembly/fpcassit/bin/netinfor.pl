#!/usr/bin/perl
use strict;

my %scaflen;
open IN, "../input/scaffold" or die "$!";
while (<IN>){
     chomp $_;
     my @unit=split("\t",$_);
     $scaflen{$unit[0]}=$unit[1];
}
close IN;
open OUT, ">../output/netinfor" or die "$!";
print OUT "Chr\tChr Length\tStart On Chr\tEnd On Chr\tHit Length Of Chr\tDistance to Last Hit\tStrand\tScaffold\tScaf Length\tStart On Scaf\tEnd On Scaf\tHit Length Of Scaf\n";
my @file=<../input/net/*>;
foreach (sort @file){
    print "$_\n";
    my $chr;
    my $chrlen;
    my $lastend=0;
    open IN, "$_" or die "$!";
        while (<IN>){
             chomp $_;
             if ($_=~/^net/){
                my @unit=split(" ",$_);
                $chr=$unit[1];
                $chrlen=$unit[2];
             } 
             if ($_=~/fill.*top/){
                 my @unit=split(" ",$_);  
                 my $refs=$unit[1];
                 my $refe=$unit[1]+$unit[2];
                 my $distance=$refs-$lastend;
                 $lastend=$refe;
                 my $qrys=$unit[5];
                 my $qrye=$unit[5]+$unit[6];
                 print OUT "$chr\t$chrlen\t$refs\t$refe\t$unit[2]\t$distance\t$unit[4]\t$unit[3]\t$scaflen{$unit[3]}\t$qrys\t$qrye\t$unit[6]\n"; 
             }
        }  
    close IN;
}
close OUT;
