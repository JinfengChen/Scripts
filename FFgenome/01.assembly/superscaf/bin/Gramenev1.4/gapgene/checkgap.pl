#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","help");


my $help=<<USAGE;
perl $0 --fasta 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $refseq=getfastaseq($opt{fasta});
foreach my $id ( sort keys %$refseq){
   while ($refseq->{$id}=~/(N+)/g){
      my $len=length $1;
      print "$id\t$len\n";
   }
}


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
 
