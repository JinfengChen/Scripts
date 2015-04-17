#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;

my %ID;

my $change_id="unknown";
$ID{$change_id}=$change_id;
my $number=1;

open (IN,$file) || die "$!";
while(<IN>){
	chomp;
	my @c=split/\t/,$_;
	if(@c > 4){
		if($c[2]=~/mRNA/){
			my $id=$1 if($c[8]=~/^ID=(\S+?);/);
			if(exists $ID{$id}){
				$change_id="$id\_$number";
				$c[8]=~s/$id/$change_id/g;
				$number++;
			}else{
				$change_id=$id;
				$ID{$id}=$id;
			}
			print join("\t",@c)."\n";
		}elsif($c[2]=~/exon/){
			my $id=$1 if($c[8]=~/^Parent=(\S+?);/);
			$c[8]=~s/$id/$change_id/g;
			print join("\t",@c)."\n";
		}
	}
}
close IN;
	
		
		
		
			
