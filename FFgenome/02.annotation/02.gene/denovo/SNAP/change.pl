#!/usr/bin/perl
use strict;
use warnings;

my $genome=shift;

open (IN,$genome) || die "$!";
$/=">";<IN>;$/="\n";
while(<IN>){
	my $scaf_name=$1 if($_=~/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp($seq);
	print ">$scaf_name\n$seq";
	$/="\n";
}
close IN;
