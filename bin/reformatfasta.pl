#!/usr/bin/perl

use FindBin qw ($Bin);
use Getopt::Long;


GetOptions(\%opt,"fasta:s","output:s","help:s");

my $help=<<USAGE;

Reform fasta format into one sequence one line, which is needed by abyss-samtoafg
USAGE

if (exists $opt{help} or !-f $opt{fasta}){
   print $help;
   exit;
}

$opt{output} ||= $opt{fasta}.".reform";


### read fasta file and output fasta if id is found in list file
$/=">";
my $id=0;
open IN,"$opt{fasta}" or die "$!";
open OUT, ">$opt{output}" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    $id++;
    my @unit=split("\n",$_);
    my $head=shift @unit;
    $head = "Scaffold".$id;
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\t//g;
    $seq=~s/\s//g;
    print OUT ">$head\n$seq\n";
 
}
close OUT;
close IN;
$/="\n";


