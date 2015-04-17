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
			$ID{$title}[$i][0]=~s/(\S+)/$1\-D2/ if($i==1);
			$ID{$title}[$i][0]=~s/(\S+)/$1\-D3/ if($i==2);
			print join("\t",@{$ID{$title}[$i]})."\n";
		}
	}
}
 
