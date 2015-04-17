#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"id:s","type:s","help");


my $help=<<USAGE;
perl $0 --id --type

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 


open IN, "$opt{id}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    print "$unit[0]\n" if ($unit[1] eq $opt{type})
}
close IN;
