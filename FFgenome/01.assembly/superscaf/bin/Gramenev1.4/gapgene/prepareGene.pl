#!/usr/bin/perl

use FindBin qw ($Bin);
use Getopt::Long;


GetOptions(\%opt,"list:s","fasta:s","output:s","help:s");

my $help=<<USAGE;
Get fasta sequences for ids in a list file
perl getidseq.pl -l id -f RAP3.fa -o nohit.fa > log  

USAGE

if (exists $opt{help} or !-f $opt{fasta}){
   print $help;
   exit;
}


### store id into hash %id
my %id;
open IN, "$opt{list}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    $id{$unit[1]}=1;    
}
close IN;


### read fasta file and output fasta if id is found in list file
$/=">";
open IN,"$opt{fasta}" or die "$!";
open OUT, ">$opt{output}" or die "$!";
open OUT1, ">$opt{output}.full" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp,2);
    my $head=$temp1[0];
    my $anno=$temp1[1];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    if (exists $id{$head}){
         #print "$head\n";
         my @temp=split(/X+/,$seq);
         #print @temp;
         my $seq1=$temp[0];
         my $seq2=$temp[$#temp];
         my $head1=$head."_5";
         my $head2=$head."_3";
         print OUT ">$head1\n$seq1\n" if ($seq1=~/\w+/);
         print OUT ">$head2\n$seq2\n" if ($seq2=~/\w+/);
         print OUT1 ">$head\n$seq\n";
    
    }
}
close OUT;
close OUT1;
close IN;
$/="\n";


