#!/usr/bin/perl

use FindBin qw ($Bin);
use Getopt::Long;


GetOptions(\%opt,"list:s","gff:s","output:s","help:s");

my $help=<<USAGE;
Get gff for ids in a list file
perl getidseq.pl -l id -g gene.gff -o id.gff > log  

USAGE

if (exists $opt{help} or keys %opt < 1){
   print $help;
   exit;
}


my $refgff=parseGFF($opt{gff});

### store id into hash %id
my %id;
open OUT, ">$opt{output}" or die "$!";
open IN, "$opt{list}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $id=$unit[0];
    #my $id=$1 if ($unit[0]=~/(.*)\_\w+/);
    if (exists $refgff->{$id}){
       print OUT "$refgff->{$id}";
    }  
}
close IN;
close OUT;



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
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return \%hash;
}


