#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","clean:s","help");


my $help=<<USAGE;
clean gene id for fasta file.
perl $0 --fasta --clean
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

cleanfasta($opt{fasta},$opt{clean});

sub cleanfasta
{
$/=">";
my %hash;
my ($file,$out)=@_;
open OUT, ">$out" or die "$!";
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    $head=$1 if ($head=~/LOC_(\w+)\.\d+/ or $head=~/(\w+)\.\d+/);
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    print OUT ">$head\n$seq\n";
}
close OUT;
close IN;
$/="\n";
}
 
