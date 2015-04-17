#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"coutdir:s","type:s","bed:s","help");


my $help=<<USAGE;
convert cout file into bed format
--coutdir: dir of cout files
--bed: file name want to output bed
perl $0 -coutdir ./ -bed FF.methylation.bed
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my @file=glob("$opt{coutdir}/*.cout");
open OUT, ">$opt{bed}" or die "$!";
foreach (@file){
   print "$_\n";
   my $file=$_;
   unless ($file=~/control/i){
      open IN, "$file" or die "$!";
      while(<IN>){
           chomp $_;
           next if ($_ eq "");
           my @unit=split("\t");
           next if ($unit[6] == 0 and $unit[7] == 0);
           next if ($unit[6] < $unit[8] or $unit[6] < 2);
           my $start=$unit[1]-1;
           my $end=$start;
           print OUT "$unit[0]\t$start\t$end\t$unit[6]\t$unit[8]\t$unit[2]\n";
      }
      close IN;
   }
}
close OUT;
