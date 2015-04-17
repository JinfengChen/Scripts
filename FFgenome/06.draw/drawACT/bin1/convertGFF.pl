#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","project:s","help");


my $help=<<USAGE;

Convert gff format file to mcscan gff format.
Example mcscan gff:
Os01 Os01g0100100 1903 9816
Sb01 Sb01g001001  2000 3000

Run: perl convertGFF.pl -gff Rice.gff -project Os
-gff: GFF file 
-project: prefix for chromosome name in mcscan gff, Os01

USAGE
 
if (keys %opt < 1){
print "$help\n";
exit();
}

open IN, "$opt{gff}" or die "$!";
open OUT, ">$opt{project}.gff" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if ($unit[2] eq "mRNA"){
         my $start=$unit[3]; 
         my $end  =$unit[4];
         $unit[0]=~s/[a-zA-Z]//g;
         my $chr =$opt{project}.$unit[0];
         my $gene;
         if ($unit[8] =~/ID=(.*);/){
             $gene=$1;
         }
         print OUT "$chr\t$gene\t$start\t$end\n";
    }
}
close IN;
close OUT;
