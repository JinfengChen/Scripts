#!/usr/bin/perl
use strict;
use warnings;

my $genome=shift;
my $gene_infor=shift;

my %scafold_name;

open (IN,$genome) || die "$!";
$/=">";<IN>;$/="\n";
while(<IN>){
	my $name=$1 if($_=~/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp($seq);
	$seq=~s/\s+//g;
	my $length=length($seq);
	if($length >= 2000){
		$scafold_name{$name}=1;
	}
	$/="\n";
}
close IN;

open (IN1,$gene_infor) || die "$!";
while(<IN1>){
	chomp;
	my @c=split;
	my $scaf_name=$c[0];
	if(exists $scafold_name{$scaf_name}){
		print join("\t",@c)."\n";
	}
}
close IN1;
