#!/usr/bin/perl

open IN, "../input/table" or die "$!";
while (<IN>){
    next if ($_ eq "");
    my @unit=split("\t",$_);
    `perl hitbes.pl $unit[0] $unit[1] > log`;
}
close IN;
