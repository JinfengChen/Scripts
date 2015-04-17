#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"coutdir:s","type:s","bed:s","help");


my $help=<<USAGE;
convert cout file into bed format, record both converted C(T) and C.
--coutdir: dir of cout files
--type: CG,CHG,CHH
--bed: dir for file want to output
perl $0 -coutdir ./ -type CG -bed FF_methylation_bed
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 
`mkdir $opt{bed}` unless (-f $opt{bed});

my @file=glob("$opt{coutdir}/*.cout");
foreach (@file){
 my $file=$_;
 unless ($file=~/control/i){
   my $scaf;
   if ($file=~/(\w+)\.cout/){
      $scaf=$1;
   }
   my $bed=$scaf.".".$opt{type}."."."bed";
   open OUT, ">$opt{bed}/$bed" or die "$!";
   print "$_\n";
   open IN, "$file" or die "$!";
      while(<IN>){
           chomp $_;
           next if ($_ eq "");
           my @unit=split("\t");
           next if ($unit[3] ne $opt{type});
           next if ($unit[6] == 0 and $unit[7] == 0);
           my $start=$unit[1]-1;
           my $end=$start;
           print OUT "$unit[0]\t$start\t$end\t$unit[6]:$unit[7]\t$unit[8]\t$unit[2]\n";
      }
   close IN;
   close OUT;
 }
}
