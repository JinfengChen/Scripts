#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blast:s","help");


my $help=<<USAGE;
perl $0 --blast

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable($opt{blast});

sub readtable
{
my ($file)=@_;
my $head=$1 if ($file=~/(.*)\.blast/);
open OUT, ">$head.pair.blast" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=~tr/\-/\_/;
    $unit[1]=~tr/\-/\_/;
    print OUT "$unit[0]\t$unit[1]\t$unit[2]\n";
}
close IN;
close OUT;
}
 
