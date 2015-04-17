#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Draw tree for all fa file in this directory
perl $0
USAGE


if ($opt{help}){
    print "$help\n";
    exit();
} 

my @file=glob("*.fa");
foreach my $fa (@file){
   next unless ($fa=~/.*\.fa$/);
   `perl /home/jfchen/FFproject/FFgenome/03.evolution/Family_analysis/bin/DrawTree.pl --protein $fa --align --nj > log 2> log2`;
}
