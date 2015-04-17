#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"masked:s","fasta:s","outfile:s","help");


my $help=<<USAGE;
perl $0 --masked --fasta --outfile  

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $seq=getfastaseq($opt{fasta});
my $mask=maskseq($opt{masked});
#####
open OUT, ">$opt{outfile}" or die "$!";
foreach my $t (keys %$seq){
    #print "$t\t$mask->{$t}->[0]\n";
    if ($mask->{$t}->[0] < 0.3){
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


sub maskseq
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
    my $length=length $seq;
    my $x=Xrate($seq);
    $hash{$head}=[$x,$length];
}
$/="\n";
return \%hash;
}
 
sub Xrate{
my ($seq)=@_;
my $length=length $seq;
my $x=$seq=~tr/atcg/atcg/;
my $Xrate=$x/$length;
return $Xrate;
}


