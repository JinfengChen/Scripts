#!/usr/bin/perl
use strict;
##remove the redundant lines

die "$0 <solar_file>\n" if(@ARGV != 1);

my $solar_file = shift;

my %uniq;
open IN,$solar_file || die "fail";
while (<IN>) {
	my @t = split /\t/;
	my $query_id = $t[0];
	my $target_id = $t[5];
	if(exists $uniq{$target_id}{$query_id}){
		delete $uniq{$target_id}{$query_id};
		next;
	}
	
	print $_;
	$uniq{$query_id}{$target_id} = 1;
}
close IN;


