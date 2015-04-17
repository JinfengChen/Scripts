#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
clean gff gene id and convert to simple format
out format:
Bradi1g02850    Bd1     1937097
perl $0 --gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

cleangff($opt{gff});

sub cleangff
{
my ($gff)=@_;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        my $head;
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*?)$/){
            $head=$1;
        }
        $head=$1 if ($head=~/LOC_(\w+)\.\d+/ or $head=~/(\w+)\.\d+/);
        print "$head\t$unit[0]\t$unit[3]\n";
    }
    

}
close IN;
}
 
