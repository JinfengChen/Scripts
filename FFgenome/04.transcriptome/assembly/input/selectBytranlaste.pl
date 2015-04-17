#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","len:s","outfile:s","help");


my $help=<<USAGE;
perl $0 --fasta --outfile
select transcript that can produce a tranlation longer by "len"
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{len} ||="100";
$opt{outfile} ||="longORF.out";

my $seq=getfastaseq($opt{fasta});
`translate $opt{fasta} > temp.translate`;
my %hash;
open IN, "temp.translate" or die "$!";
while (<IN>){
    my $line=$_;
    chomp $line;
    if ($line=~/^> (.*)\.\d+\s+length (\d+)\,/){
       #print "$1\t$2\n";
       $hash{$1}= $hash{$1} < $2 ? $2 : $hash{$1};
    }
}
close IN;

open OUT, ">$opt{outfile}" or die "$!";
foreach my $t (keys %hash){
    #print "$t\t$hash{$t}\n";
    if ($hash{$t} >= $opt{len}){
       print OUT ">$t\n$seq->{$t}\n";
    }
}


####
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
    #print "$head\n";
    $hash{$head}=$seq;
}
close IN;
$/="\n";
return \%hash;
}

