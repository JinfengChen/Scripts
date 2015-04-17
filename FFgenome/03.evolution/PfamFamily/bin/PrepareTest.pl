#!/usr/bin/perl
use Getopt::Long;
use warnings;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Read Pfamfamily.txt, generate a file to run Chi Square Test using ChiSqTest.pl.
Input File:
Pfam    OB      OS      Total   Annotation
PF00004 103     106     209     ATPase, AAA-type, core
PF00005 119     110     229     ABC transporter-like

Output File:
Pfam    OB      OB_other        OS      OS_other        Total   Annotation
PF00004 103     23738   106     24915   209     ATPase, AAA-type, core
PF00005 119     23722   110     24911   229     ABC transporter-like
PF00006 8       23833   6       25015   14      ATPase, F1/V1/A1 complex, alpha/beta subunit, nucleotide-binding domain
PF00008 0       23841   1       25020   1       EGF-like


Run: perl PrepareTest.pl Pfamfamily.txt Pfamfamily4test.txt
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
}


open IN, "$ARGV[0]" or die "$!";
######determine the number of species########
my $head=<IN>;
chomp $head;
my @head=split("\t",$head);
my $last;
for(my $i=0;$i<@head;$i++){
   if ($head[$i] eq "Total"){
        $last=$i-1;
   }
}
#############################################
my @total; ## total for each species column
my @record;## record for each Pfam
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   push (@record,[@unit]);
   for(my $i=1;$i<=$last;$i++){
      $total[$i-1]+=$unit[$i];
      #print "$i\t$total[$i-1]\n";
   }
}
close IN; 

open OUT, ">$ARGV[1]" or die "$!";
my $newh=$head[0];
for(my $i=1;$i<=$last;$i++){
   my $other=$head[$i]."_"."other";
   $newh.="\t$head[$i]\t$other";
}
$newh.="\t$head[$last+1]\t$head[$last+2]";
print OUT "$newh\n";
for(my $j=0;$j<@record;$j++){
    print OUT "$record[$j][0]\t";
    for(my $i=1;$i<=$last;$i++){
       my $left=$total[$i-1]-$record[$j][$i];
       print OUT "$record[$j][$i]\t$left\t";
    }
    print OUT "$record[$j][$last+1]\t$record[$j][$last+2]\n";
}
close OUT;
