#!/usr/bin/perl
#### clear tigr rice gff3 file to remove repeat related gene and alternative splicing gene

my %temodel;
open IN, "../data/gff/all.TE_related" or die "$!";
while (<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\t",$_);
     unless (exists $temodel{$unit[0]}){
         $temodel{$unit[0]}++;
     }
}
close IN;
my $number=keys %temodel;

$/="###";
my $counter;
open OUT, ">../data/gff/all.con.gff3.clear" or die "$!";
open IN, "../data/gff/all.con.gff3" or die "$!";
while(<IN>){
     chomp $_;
     next if (length $_ < 5);
     if ($_=~/MSU_osa1r6	gene\s+.*Alias\=(LOC_Os\d{2}g\d{5})\.\d{1}\n/s){
          #print "$1\n";
          unless (exists $temodel{$1}){
               $counter++;
               print OUT "$_###";
          }
     }
}
close IN;
close OUT;

print "Final Gene Model: $counter\n";
