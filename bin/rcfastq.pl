#!/usr/bin/perl

# AUTHOR: Joseph Fass
# LAST REVISED: September 2010
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2008 The Regents of University of California, Davis Campus.
# All rights reserved.

# USAGE: cat sequence.fastq | ./rcFastq.pl > sequence.rc.fastq

while (<>) {
	$h1 = $_;
	$s = <>;
	$h2 = <>;
	$q = <>;
	chomp $s;
	chomp $q;
	$s = reverse $s;
	$s =~ tr/ATCGNatcgn/TAGCNtagcn/;
	$q = reverse $q;
	print $h1.$s."\n".$h2.$q."\n";
} # while


