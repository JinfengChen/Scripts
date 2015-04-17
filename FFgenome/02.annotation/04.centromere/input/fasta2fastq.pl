#!/usr/bin/perl
use strict;
use warnings;
die ("Usage: fasta2fastq <fasta.file> <qual.file>") unless  (scalar @ARGV) == 2;
open FASTA, $ARGV[0] or die "cannot open fasta: $!\n";
open QUAL, $ARGV[1] or die "cannot open qual: $!\n";
my $offset = 33; # I think this was 33 for sanger FASTQ, change this if required!
my $count = 0;
local($/) = "\n>"; # split the input fasta file by FASTA records
# this is some splitting of the fasta by line
while (my $fastarec = <FASTA>) {  
   chomp $fastarec;
   my ($fid, @seq) = split "\n", $fastarec;
   my $seq = join "", @seq; $seq =~ s/\s//g;
   my $qualrec = <QUAL>;
   chomp $qualrec;
   my ($qid, @qual) = split "\n", $qualrec;
   @qual = split /\s+/, (join( " ", @qual));  # convert score to character code:
   my @qual2 = map {chr($_+$offset)} @qual;
   my $quals = join "", @qual2;
   die "missmatch of fasta and qual: '$fid' ne '$qid'" if $fid ne $qid;
   $fid =~ s/^\>//;
   print STDOUT (join( "\n", "@".$fid, $seq, "+$fid", $quals), "\n");
   $count++;
}
close (FASTA);
close (QUAL);
print STDERR "wrote $count entries\n";

