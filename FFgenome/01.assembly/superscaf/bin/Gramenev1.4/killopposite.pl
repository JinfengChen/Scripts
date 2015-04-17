#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"run:s","help");


my $help=<<USAGE;
perl $0 --run 1

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my %kill;
open IN, "opposite.list" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   $kill{$unit[0]}=1;
}
close IN;

open OUT, ">Gramene.final.v1.4.pep.iprscan" or die "$!";
open IN, "OBa.gramene.final.v1.4.pep.iprscan" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   unless (exists $kill{$unit[0]}){
      print OUT "$_\n";
   }else{
      $count++;
      print "$_\n";
   }
}
close IN;
close OUT;

print "$count\n";
