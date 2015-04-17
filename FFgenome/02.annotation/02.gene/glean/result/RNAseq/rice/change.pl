#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;

my %gene;

open (IN,$file) || die "$!";
while(<IN>){
	my @c=split/\t/,$_;
	if($c[2]=~/mRNA/){
		my $id=$1 if($c[8]=~/^ID=(\S+?);/);
		push @{$gene{$id}{mRNA}},(@c);
	}
	if($c[2]=~/CDS/ || $c[2]=~/exon/){
		my $id=$1 if($c[8]=~/^Parent=(\S+)/);
		push @{$gene{$id}{CDS}},[@c];
	}
}
close IN;

foreach my $title(sort keys %gene){
	print "$gene{$title}{mRNA}[0]\t$gene{$title}{mRNA}[1]\tmRNA\t$gene{$title}{mRNA}[3]\t$gene{$title}{mRNA}[4]\t$gene{$title}{mRNA}[5]\t$gene{$title}{mRNA}[6]\t$gene{$title}{mRNA}[7]\tID=$title;\n"; 
	my $number=@{$gene{$title}{CDS}};
	for(my $i=0;$i<$number;$i++){
		print "$gene{$title}{CDS}[$i][0]\t$gene{$title}{CDS}[$i][1]\tCDS\t$gene{$title}{CDS}[$i][3]\t$gene{$title}{CDS}[$i][4]\t$gene{$title}{CDS}[$i][5]\t$gene{$title}{CDS}[$i][6]\t$gene{$title}{CDS}[$i][7]\tParent=$title;\n";
	}
}
		
