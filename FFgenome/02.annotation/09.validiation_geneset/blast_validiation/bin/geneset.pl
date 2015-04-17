#!/usr/bin/perl

use FindBin qw ($Bin);

my $input="$Bin/../input";
my $output="$Bin/../output";

#####store annotation infor of RAP gene in hash %note
my %note;
open IN, "$input/RAP3inf" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    $note{$unit[0]}="$unit[1]\t$unit[6]";
}
close IN;

#####parse blast tab file, output gene infor for which has no hit in FF gene
my %hit;
open IN, "$output/RAP3nrORF.fa.blast.tab" or die "$!";
while (<IN>){
      chomp $_;
      next if ($_ eq "");
      next if ($_ =~/Query_id/);
      my @unit=split("\t",$_);
      my $hitp=abs ($unit[3]-$unit[2])/$unit[1]; 
      if ($unit[8] >= 0.5){
          $hit{$unit[0]}=1;   
      }
}
close IN;

open OUT, ">$output/RAPnohit" or die "$!";
foreach (sort keys %note){
        unless (exists $hit{$_}){
           print OUT "$_\t$note{$_}\n";    
        }
}
close OUT;



