#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blastm8:s","project:s","help");

my $help=<<USAGE;
Convert Blastm8 file from ortholog_paralog_pipeline to cynetenator blast file.

Example mcscan blast file (gene1 gene2 p-value); 
Ab003491        Os08g0135900    score
Ab013449        Os08g0539400    score
Ab013449        Os08g0539700    score
Ab013449        Os12g0559500    score
Ab013449        Os09g0517200    score


Run: perl $0 -blastm8 all_vs_all.blast.m8 -project 3way
-blast: blasttable file
-project: prefix for result blast file

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}


my $python="python";
my $filter_blast="/home/jfchen/FFproject/tools/mcscan0.8/filter_blast.py";

my $unfiltered=$opt{blastm8}.".unfiltered";
my $blast =$opt{project}.".blast";
#system "cut -f 1,2,11 $opt{blast} > $unfiltered";
& name($opt{blastm8},$unfiltered);
##filter self matches and unique multi hsp
##
system "$python $filter_blast $unfiltered $blast";


sub name
{
my ($file,$out)=@_;
open IN, "$file" or die "$!";
open OUT, ">$out" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   $unit[0]=$1 if ($unit[0]=~/(.*)\_\w+$/);
   $unit[1]=$1 if ($unit[1]=~/(.*)\_\w+$/);
   print OUT "$unit[0]\t$unit[1]\t$unit[11]\n";
}
close IN;
close OUT;
}


