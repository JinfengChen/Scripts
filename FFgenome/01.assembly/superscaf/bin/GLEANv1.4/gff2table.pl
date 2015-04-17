#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,,"gff:s","help");


my $help=<<USAGE;

Run: perl $0 -gff rice.gene.gff > rice.table 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

parseGFF($opt{gff});

###########################################

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        print "$unit[0]\t$id\t$unit[3]\t$unit[4]\n";
    }

}
close IN;
}

