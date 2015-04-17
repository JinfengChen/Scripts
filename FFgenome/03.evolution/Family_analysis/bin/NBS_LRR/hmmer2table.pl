#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"hmmer:s","help");


my $help=<<USAGE;
Parse hmmer2.3.2 result of hmmpfam into table.
perl $0 --hmmer test.hmmsearch > test.hmmsearch.table
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

parsehmmer($opt{hmmer});

sub parsehmmer
{
my ($file)=@_;
#print "$file\n";
$/="\/\/\n";
open IN, "$file" or die "$1";
while(<IN>){
    my $record=$_;
    #print "AAA:$record\n";
    if ($record=~/Query sequence\: (.*?)\n[\d\D]+\s---\n(.*)\n([\d\D]+)\s+-------\n([\d\D]+)\n\nAlignments of top-scoring domains/){
       my $gene=$1;
       my @hit=split("\n",$4);
       foreach my $h (@hit){
          my @domain=split(" ",$h);
          print "$gene\t$domain[0]\t$domain[2]\t$domain[3]\t$domain[5]\t$domain[6]\t$domain[8]\t$domain[9]\n";
       }
    }
}
close IN;
$/="\n";
}
