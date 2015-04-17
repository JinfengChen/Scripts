#!/usr/bin/perl
use strict;
use FindBin qw ($Bin);
## parse blastout and construct the relationship between scaffold and physical contig
## read in bes length and unmasked repeat length stored in hash
 if (@ARGV<1){
       print "Usage: perl scaf2ctg.pl blastoutfilem8 > log &\n";
       print "This scripts need bes.length file and clone, which contain infomration of bes:\n";
       print "bes.length:\nClone\tlength\tUnmasked Length\nOB__Ba0001A01.f	603	211\n";
       print "clone:\nClone\tContig\tStart CB\tEnd CB\na0001A01	ctg309	248	348\n";
       print "The result file is sort.scafmap, which is sorted based scaf and then for ctg and for position\n";
       print "Scaffold\tClone\tContig\tContig Positon (CB unit)\tScaffold Position (bp)\n";
       print "chr04.con10360001:11694001	a0038I13	121	126	723604.5\n";
       exit;
 }
my $input="$Bin/../input";
my $blastfile=$ARGV[0];
my %beslength;
my %repeat;
open BES, "$input/bes.length" or die "can not open my bes";
while (<BES>){
      chomp $_;
      my @unit=split("\t",$_);
      my $indicator=$unit[2]/$unit[1];
      if ($indicator > 0.6){
          $repeat{$unit[0]}=0; 
      }else{
          $repeat{$unit[0]}=1;
      }
      $beslength{$unit[0]}=$unit[1];
}
close BES;
## read in clone inf store contig and position of each clone on contig in hash
my %clone2ctg;
my %clone2pos;
my $multictg;
open CLONE, "$input/clone" or die "can not open my clone file";
    while (<CLONE>){
          chomp $_;
          my @unit=split("\t",$_);
          unless (exists $clone2ctg{$unit[0]}){
                 $clone2ctg{$unit[0]}=$unit[1];
                 $clone2pos{$unit[0]}=$unit[2]+($unit[3]-$unit[2])/2;
          }else{
                 $multictg++;
          }
    }
    #print "$multictg\n"; 
close CLONE;

## read blast out file, get the best hit for each clone.
## use identify and coverage as cutoff
## write the result to file scafmap, scaffold\tclone\tcontig\tposition\n, and sort by scaffold then by contig and by position.
my %besthit;
my %hitonscafpos;
open BLAST, "$blastfile" or die "can not open my blast file";
       while (<BLAST>){
            chomp $_;
            my @unit=split("\t",$_);
            $unit[0]=~/OB__B(\w+)\.\w{1}/;
            my $clone=$1;
            my $hitlength=abs($unit[7]-$unit[6]);
            my $hitonscaf=$unit[8]+($unit[9]-$unit[8])/2;
            my $coverage=$hitlength/$beslength{$unit[0]}; ### coverage on bes 
            unless ($besthit{$clone} or $unit[2] < 95 or $coverage < 0.6 or $repeat{$unit[0]}==1){ ## identity < 85, repeat
                   $besthit{$clone}=$unit[1];
                   $hitonscafpos{$clone}=$hitonscaf;
            }
       }
close BLAST;
open OUT, ">scafmap" or die "can not open my scafmap";
foreach (sort {$besthit{$a} cmp $besthit{$b}} keys %besthit){ # sort the hash by value 
       print OUT  "$besthit{$_}\t$_\t$clone2ctg{$_}\t$clone2pos{$_}\t$hitonscafpos{$_}\n"; 

}
close OUT;
system "msort -k '1,n3,n4' scafmap > sort.scafmap &";
print "Done!";




