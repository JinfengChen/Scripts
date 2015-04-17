#!/usr/bin/perl 
use strict;
use warnings;

my $file1=shift;
my $file2=shift;

my %ID;
my %seq;

open (IN1,$file1) || die "$!";
while(<IN1>){
	my @c=split/\t/,$_;
	if($c[2]=~/mRNA/){
		my $id=$1 if($c[8]=~/ID=(\S+);/);
		$ID{$id}=$id;
	}
}
close IN1;

open (IN2,$file2) || die "$!";
$/=">";<IN2>;$/="\n";
while(<IN2>){
	my $seq_name=$1 if($_=~/^(\S+)/);
	$/=">";
	my $sequence=<IN2>;
	chomp($sequence);
	$sequence=~tr/atcgu/ATCGU/;
	if(exists $ID{$seq_name}){
		print ">$seq_name\n$sequence";
	}
	$/="\n";
}
close IN2;

	
