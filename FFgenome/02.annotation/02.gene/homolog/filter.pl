#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;

my %ID;
open IN,$file || die "$!";
while(<IN>){
	chomp;
	my @c=split/\t/,$_;
	my $id=$c[0];
	push @{$ID{$id}},[@c];
}
close IN;

foreach my $title(keys %ID){
	@{$ID{$title}}=reverse(sort {$a->[10]<=>$b->[10]}@{$ID{$title}});
	for(my $i=0;$i<3;$i++){
		if(defined @{$ID{$title}[$i]}){
			print join("\t",@{$ID{$title}[$i]})."\n";
		}
	}
}
 
