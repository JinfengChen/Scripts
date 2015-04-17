#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"block:s","help");


my $help=<<USAGE;
Summary collinearity from block file
perl $0 --block
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $summary= $1 if ($opt{block}=~/(.*)\.blocks/);

& parseblock($opt{block});

sub parseblock
{
my ($block)=@_;
my ($refer,$refn,$genen);
my (%ref,%qry,$pair);
open OUT, ">$block.table" or die "$!";
open IN, "$block" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my $line=$_;
     if ($line=~/## View (\d+)\: pivot (\w+)$/){
        $refn=$1;
        $refer=$2; 
     }elsif($line=~/^#/){
        next;
     }elsif($line=~/^\s+(\d+)\-\s*(\d+)\:\s+?(.*)$/){
        my $genen=$2;
        my $genes=$3;
        my @unit=split(" ",$genes);
        #$genen=$1 if ($unit[1]=~/(\d+)\:/);
        #print "$refer\t$refn\t$genen\t$unit[0]\t$unit[1]\n" if ($unit[0]=~/\w+/ and $unit[1]=~/\w+/);
        if ($unit[0]=~/\w+/ and $unit[1]=~/\w+/){
            $pair++;
            print OUT "$refer\t$refn\t$genen\t$unit[0]\t$unit[1]\n";
            my @gene1=split(";",$unit[0]);
            my @gene2=split(";",$unit[1]);
            foreach my $g1 (@gene1){
               $ref{$g1}=1;
            }
            foreach my $g2 (@gene2){
               $qry{$g2}=2;
            }
        }
        
     } 
}
close IN;
close OUT;
my $refgenen=keys %ref;
my $qrygenen=keys %qry;
print "Total Pair: $pair\n";
print "Reference Gene Number: $refgenen\n";
print "Target Gene Number: $qrygenen\n";
}## end of sub
