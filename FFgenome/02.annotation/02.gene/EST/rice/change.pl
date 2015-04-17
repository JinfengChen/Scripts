#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $file=shift;

open (IN,$file) || die "$!";
$/=">";<IN>;$/="\n";
while(<IN>){
	my $id=$1 if($_=~/^gnl\|UG\|Os\#(\S+)/);
	 print ">$id\n";
	$/=">";
	my $seq=<IN>;
	chomp($seq);
	$seq=~tr/atcgn/ATCGN/;
	my @c=split(/\n/,$seq);
#	print Dumper @c;
#	exit;
	my $number=@c;
	for(my $i=0;$i<$number;$i++){
		next if($c[$i]=~/\#/ || !$c[$i]);
		print "$c[$i]\n";
	}
	$/="\n";
}
close IN;
	
