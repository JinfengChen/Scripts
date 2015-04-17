#!/usr/bin/perl
use Getopt::Long;
use warnings;
use strict;

my %opt;

GetOptions (\%opt,"table:s","share","help");


my $help=<<USAGE;

Run: perl sumTable.pl -table OBa.family.matrix -share > OBa.family.table
USAGE

my %name=(
     "osa" => "O. sativa",
     "oba" => "O. brachyantha",
     "brd" => "B. distachyon",
     "zma" => "Z. mays",
     "sbi" => "S. bicolor",
     "ptr" => "P. trichocarpa"
);



if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

open IN, "$opt{table}" or die "$!";
open OUT, ">temp" or die "$!";
my $first=<IN>;
chomp $first;
my @species=split("\t",$first);
my $second=<IN>;
chomp $second;
my @total=split("\t",$second);
while(<IN>){
   print OUT "$_";
}
close OUT;
close IN;

print "Species\t#genes\t#families\t#gene in families\tUnique Families\tUnique Genes\tSingle-gene Families\t#genes per Family\n";
for(my $i=1;$i<=@species;$i++){
   print "$name{$species[$i-1]}\t";
   print "$total[$i-1]\t";
   ### number of family for each species
   my $family="awk '\$$i > 0' temp | wc -l";
   my $familyn=`$family`;
   chomp $familyn;
   print "$familyn\t";
   
   ### number of gene in family for each species
   my $geneinfam="awk '{if (\$$i> 0);sum+=\$$i}END{print sum}' temp";
   my $geneinfamn=`$geneinfam`;
   chomp $geneinfamn;  
   print "$geneinfamn\t";
  
   ### number of shared family in all species
 
   ### number of uniq family for each species
   my @temp;
   for(my $j=1;$j<=@species;$j++){
      if ($j == $i){
          push (@temp,"\$$j > 0");
      }else{
          push (@temp,"\$$j == 0");
      }
   }
   my $temp=join("&&",@temp);
   my $uniqfam="awk '$temp' temp | wc -l";
   my $uniqfamn=`$uniqfam`;
   chomp $uniqfamn;
   print "$uniqfamn\t"; 


   ### number of gene in uniq family for each species
   my $uniqfile="awk '$temp' temp > uniq";
   `$uniqfile`;
   my $uniqgene="awk '{sum+=\$$i}END{print sum}' uniq"; 
   my $uniqgenen=`$uniqgene`;
   chomp $uniqgenen;
   print "$uniqgenen\t";
   
   ### number of single gene or family for each species
   my $single="awk '\$$i == 1' temp | wc -l";
   my $singlen=`$single`;
   chomp $singlen;
   print "$singlen\t";

   ### number of gene per family
   my $geneperfamily=sprintf("%.2f",$geneinfamn/$familyn);
   print "$geneperfamily\n";
   
   if ($opt{share}){ ### use this when there are four species
   ### share between three species
   my @temp3;
   for(my $j=1;$j<=@species;$j++){
      if ($j == $i){
          push (@temp3,"\$$j == 0");
      }else{
          push (@temp3,"\$$j > 0");
      }
   }
   my $temp3=join(" && ",@temp3);
   my $share3="awk '$temp3' temp | wc -l";
   my $share3n=`$share3`;
   chomp $share3n;
   my $share3out="no-$species[$i-1]:$share3n";
   print "$share3out\t";
   

   ### share between two species,unfinished
   for(my $j=1;$j<=@species;$j++){
      if ($j != $i){
         my @temp2;
         for (my $a=1;$a<=@species;$a++){
            if ($a == $j or $a == $i){
               push (@temp2,"\$$a > 0");
            }else{
               push (@temp2,"\$$a == 0");
            }
         }
         my $temp2=join("&&",@temp2);
         my $share2="awk '$temp2' temp | wc -l";
         my $share2n=`$share2`;
         chomp $share2n;
         my $share2out="$species[$i-1] $species[$j-1]:$share2n";
         print "$share2out\t";
      }
   }
   print "\n";
   ### end if share
   }
}

   ### number of shared family in all species
   my @temp;
   for(my $j=1;$j<=@species;$j++){
      push (@temp,"\$$j > 0");
   }
   my $temp=join("&&",@temp);
   my $sharedfam="awk '$temp' temp | wc -l";
   my $sharedfamn=`$sharedfam`;
   chomp $sharedfamn;
   print "Shared in All: $sharedfamn\n";


