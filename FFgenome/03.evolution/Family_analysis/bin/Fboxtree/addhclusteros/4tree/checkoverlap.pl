#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dir:s","help");


my $help=<<USAGE;


USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash;
my @file=glob("$opt{dir}/*.fa");
foreach my $file (@file){
   getfastaseq($file,\%hash);
}

foreach(sort keys %hash){
   my $type=join(":",@{$hash{$_}});
   print "$_\t$type\n" if @{$hash{$_}} > 0;
}


sub getfastaseq
{
$/=">";
my ($file,$ref)=@_;
my $type=$1 if ($file=~/$opt{dir}\/(.*)\.fa$/);
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
    push (@{$ref->{$head}},$type);
}
$/="\n";
}
 
