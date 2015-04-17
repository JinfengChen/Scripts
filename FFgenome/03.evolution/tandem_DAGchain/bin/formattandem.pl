#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"tandem:s","project:s","help");


my $help=<<USAGE;
perl $0 --tandem 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$/="##";
my $count;
open OUT1, ">$opt{project}.tandem.id.txt" or die "$!";
open OUT2, ">$opt{project}.tandem.family.txt" or die "$!";
open IN, "$opt{tandem}" or die "$!";
while(<IN>){
  chomp $_;
  next if ($_=~/^$/);
  $count++;
  my @unit=split("\n",$_);
  shift @unit;
  my %hash;
  foreach my $pair (@unit){
     my @array=split("\t",$pair);
     $hash{$array[1]}=1;
     $hash{$array[5]}=1;
  }
  my @gene=keys %hash;
  my $gene=join("\t",@gene);
  print OUT2 "$gene\n";
  foreach(keys %hash){
     print OUT1 "$_\n";
  }
}
close IN;
close OUT1;
close OUT2;
$/="\n";

