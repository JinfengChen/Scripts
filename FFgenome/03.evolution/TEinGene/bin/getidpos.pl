#!/usr/bin/perl

use FindBin qw ($Bin);
use Getopt::Long;


GetOptions(\%opt,"list:s","gff:s","output:s","help:s");

my $help=<<USAGE;
Get position for ids in a list file
perl getidseq.pl -l id -g gene.gff -o nohit.fa > log  

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
    $unit[0]=$1 if ($unit[0]=~/(.*)_E\d+/);
    $id{$unit[0]}=1;
    if (exists $refgff->{$unit[0]}){
       print OUT "$unit[0]\t$refgff->{$unit[0]}->[0]\t$refgff->{$unit[0]}->[1]\t$refgff->{$unit[0]}->[2]\t$refgff->{$unit[0]}->[3]\n";
    }  
}
close IN;
close OUT;



#######################
sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $unit[0]=~s/chr0//;
        $unit[0]=~s/chr//;
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $hash{$id}=[$seq,$unit[3],$unit[4],$unit[6]];
    }
}
close IN;
return \%hash;
}


