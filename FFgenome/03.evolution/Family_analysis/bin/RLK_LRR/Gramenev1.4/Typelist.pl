#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dir:s","help");


my $help=<<USAGE;
perl $0 --dir ./
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my @file=glob("$opt{dir}/*.fa");
foreach my $file (@file){
    if ($file=~/.*\.(.*?)\.fa$/){
       my $type=$1;
       my $refid=getfastaid($file);
       foreach my $g(keys %$refid){
          print "$g\t$type\n";
       }
    }
}

sub getfastaid
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
    $hash{$head}=1;
}
$/="\n";
return \%hash;
}

