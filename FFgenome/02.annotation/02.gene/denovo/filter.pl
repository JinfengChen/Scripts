#!/usr/bin/perl
use strict;
use warnings;

my $file1=shift;
my $file2=shift;

my %ID;
open (IN1,$file1) || die "$!";
while(<IN1>){
	chomp;
	my @c=split;
	my $id=$c[0];
	$ID{$id}=$id;
}
close IN1;


open (IN2,$file2) || die "$!";
while(<IN2>){
	chomp;
	my @c=split/\t/,$_;
	if($c[8]=~/ID=(\S+?);/ || $c[8]=~/Parent=(\S+?);/){
		my $id=$1;
		if(not exists $ID{$id}){
			print join("\t",@c)."\n";
		}
	}
}

			
