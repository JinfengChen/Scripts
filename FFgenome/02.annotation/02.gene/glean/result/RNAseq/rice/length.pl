#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;
open (IN,$file) || die "$!";
$/=">";<IN>;$/="\n";
while(<IN>){
	my $id=$1 if($_=~/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp($seq);
	$seq=~s/\s+//g;
	my $length=length($seq);
	print "$id\t$length\n";
	$/="\n";
}
close IN;
