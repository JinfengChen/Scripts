#!/usr/bin/perl

use FindBin qw ($Bin);
use Getopt::Long;


GetOptions(\%opt,"list:s","iprscan:s","output:s","help:s");

my $help=<<USAGE;
Get iprscan for ids in a list file
perl $0 -l id -i gene.iprscan -o id.iprscan > log  

USAGE

if (exists $opt{help} or keys %opt < 1){
   print $help;
   exit;
}


my $ref=readtable($opt{list});

open IN, "$opt{iprscan}" or die "$!";
open OUT, ">$opt{output}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if (exists $ref->{$unit[0]}){
       print OUT "$_\n";
    }
}
close IN;
close OUT;

#################

sub readtable
{
### store id into hash %id
my %hash;
open IN, "$opt{list}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_); 
    my $id;
    if ($unit[0]=~/(.*)\_\w+$/){
       $id=$1;
    }else{
       $id=$unit[0];
    }
    $hash{$id}=1;
    #print "$id\n";
}
close IN;
return \%hash;
}
