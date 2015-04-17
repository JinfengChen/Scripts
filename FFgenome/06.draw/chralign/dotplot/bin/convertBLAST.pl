#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blast:s","project:s","help");

my $help=<<USAGE;
Convert Blasttable to mcscan blast file.

Example mcscan blast file (gene1 gene2 p-value); 
Ab003491        Os08g0135900    0.0
Ab013449        Os08g0539400    0.0
Ab013449        Os08g0539700    0.0
Ab013449        Os12g0559500    4e-130
Ab013449        Os09g0517200    8e-123


Run: perl convertBLAST.pl -blast ob_os.blasttable -project ob_os 
-blast: blasttable file
-project: prefix for result blast file

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}


my $python="python";
my $filter_blast="/home/jfchen/FFproject/tools/mcscan0.8/filter_blast.py";

my $unfiltered=$opt{blast}.".unfiltered";
my $blast =$opt{project}.".blast";
system "cut -f 1,5,14 $opt{blast} > $unfiltered";
system "$python $filter_blast $unfiltered $blast";


