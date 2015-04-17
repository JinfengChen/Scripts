#!/usr/bin/perl
use strict;
use warnings;

my $file1=shift;
my $file2=shift;

my %scaff;

open (IN1,$file1) || die "$!";
while(<IN1>){
    my @c=split;
    my $scaf_name=$c[0];
    $scaff{$scaf_name}=$scaf_name;
}
close IN1;

my %seq;

open (IN2,$file2) || die "$!";
$/=">";<IN2>;$/="\n";
while(<IN2>){
    chomp;
    my $scaf_name=$1 if($_=~/^(\S+)/);
    $/=">";
    my $sequence=<IN2>;
    chomp($sequence);
    #$sequence=~s/\s+//g;
    $seq{$scaf_name}=$sequence;
    $/="\n";
}
close IN2;

foreach my $title(sort keys %seq){
    if(exists $scaff{$title}){
        print ">$title\n$seq{$title}";
    }
}
    
