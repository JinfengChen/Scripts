#!/usr/bin/perl
#
# edit $/ = ">";
# read fasta and get length
#

use strict;
use warnings;

# help
die "Usage: perl $0 <fasta>\n" if (@ARGV != 1);
my (@len,$total_len);
open(FA,$ARGV[0]) or die "fail to read $ARGV[0]\n";
$/ = ">";
while(<FA>){
	chomp;
	next if ($_ eq "");
	my $name = (split(/\s/,$_))[0];
	my $seq = (split(/\n/,$_,2))[1];
	$seq =~ s/\n//g;
	$seq =~ s/\s//g;
	$seq =~ s/\r//g;
	#print ">$name 51-100\n" . substr($seq,50,50) . "\n";
	push(@len,length($seq));
	$total_len += length($seq);
	#print "$name\t" . length($seq) . "\n";
}
close FA;

#print "length list: @len\n";

my $sum;
my $len_90 = $total_len*0.9;
#print "*** $len_50\n";
foreach (sort{$b <=> $a} @len){
	#print "** $_\n";
	$sum += $_;
	if($sum >= $len_90){
		print "$_";
                #print "N50: $_\n";
		exit();
	}
}


