#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"orth:s","seq1:s","seq2:s","help");


my $help=<<USAGE;
Read *.orth file, get fasta sequence for each pair into a fa file.
The output will be in a directory named as *.

Run: perl PrepareOrthFas.pl -orth OB-OS.orth -seq1 OB -seq2 OS
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $refseq1=getfastaseq($opt{seq1});
my $refseq2=getfastaseq($opt{seq2});

my $bn=`basename $opt{orth}`;
print "$bn\n";
if ($bn=~/(.*)\.orth/){
    mkdir $1 unless (-e $1);
    open IN, "$opt{orth}" or die "$!";
    while(<IN>){
        chomp $_;
        my @unit=split("\t",$_);
        my $fa=$1."/$unit[0].fa";
        open OUT, ">$fa" or die "$!";
            print OUT ">$unit[0]\n$refseq1->{$unit[0]}\n";
            print OUT ">$unit[1]\n$refseq2->{$unit[1]}\n";
        close OUT;
    }
    close IN;
}


##########sub function#################
sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    $seq=~s/\r//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}




 
