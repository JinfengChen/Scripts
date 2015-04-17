#!/usr/bin/perl
=header
Convert RPKM to mean-per-bp-coverage across gene.
perl $0 -rpkm RAP3.rpkm > log 2> log2 &
-rpkm: rpkm table
=cut

use Getopt::Long;
use warnings;
use strict;

our %opt;
GetOptions(\%opt,"rpkm:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -rpkm RAP3.rpkm\n";
   exit();
}


#OBa 30526202
#rice 16172828 
$totalread="30526202"; 

