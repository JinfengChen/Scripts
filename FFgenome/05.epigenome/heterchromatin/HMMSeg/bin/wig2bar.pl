#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;
use warnings;
#use strict;

our %opt;

GetOptions(\%opt,"wig:s","help");

my $help=<<USAGE;
Change the wig result of HMMSeg to bar file, which can be view in cisgenome.
perl wig2bar.pl -wig ../input/rice_methylation_wave
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}

my @file=glob("$opt{wig}/*.bed.wig");
foreach my $file (@file){
    print "$file\n";
    my $head;
    if ($file=~/(.*)\.bed\.wig/){
       $head=$1;
    }
    deal($file,$head);
}


######
sub deal
{
my ($file,$head)=@_;
my $tag=0;
#open OUT1, ">$head.wig.bar.txt" or die "$!";
#open OUT2, ">$head.wig.smooth.bar.txt" or die "$!";
open OUT3, ">$head.wig.segment.txt" or die "$!";
my $fh;
my @collection;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   if ($_=~/^browser/){
      next;
   }elsif($_=~/name=(data_(\d+))/){
      $tag=1;
      $fh=$1;
      #print "Data:$fh\t$2\n";
      push (@collection,$fh);
      open $fh, ">$head.wig.$2.bar.txt" or die "$!";
   }elsif($_=~/name=(smoothed_data_(\d+))/){
      $tag=2;
      $fh=$1;
      #print "Smooth_Data:$fh\t$2\n";
      push (@collection,$fh);
      open $fh, ">$head.wig.$2.smooth.bar.txt" or die "$!";
   }elsif($_=~/name=viterbi_segmentation/){
      $tag=3;
   }elsif($tag==1){
      my @unit=split("\t",$_);
      #print "TAG1:$fh\n";
      print $fh "$unit[0]\t$unit[1]\t$unit[3]\n";
   }elsif($tag==2){
      my @unit=split("\t",$_);
      #print "TAG2:$fh\n";
      print $fh "$unit[0]\t$unit[1]\t$unit[3]\n";
   }elsif($tag==3){
      #print "TAG3\n";
      print OUT3 "$_\n";
   }
}
close IN;
#close OUT1;
#close OUT2;
close OUT3;
foreach(@collection){
   close $_;
}

}


