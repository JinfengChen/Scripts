#!/usr/bin/perl

## this program is designed to find BES pairs that connecting two neighbor scaffold.
## Useage: perl hitbes.pl Scaffold000001 Scaffold000002
## the result will contain Fhit, Rhit, which are total hits number in all scaffold for f and r BES. 


use strict;

my %bes;
open IN, "../input/bes.length" or die "$!";
while (<IN>){
     chomp $_;
     my @unit=split("\t",$_);
     if ($unit[0]=~/OB__B(\w+)\.(\w+)/){
          if (exists $bes{$1}){
              $bes{$1}.="\t$2\t$unit[1]\t$unit[2]";
          }else{
              $bes{$1}="$2\t$unit[1]\t$unit[2]";
          }
     }
}
close IN;


my %ref;
my %qry;
my %ctg;
my %beshit;
my $flag;
open IN, "../input/FFversion2_BSS_results_super-scaffold.brachyanthaBES.filter.bss" or die "$!";
while (<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    if ($_=~/#Hit table/){
       <IN>;
       $flag=1;
    }elsif ($flag){
      if (exists $beshit{$unit[3]}){
        my $rev=substr($unit[0],8,1);
        #print "$rev\n";
        if ($rev eq "f"){
           my $f=$beshit{$unit[3]}->[0];
           $f++;
           $beshit{$unit[3]}->[0]=$f;
        }else{
           my $r=$beshit{$unit[3]}->[1];
           $r++;
           $beshit{$unit[3]}->[1]=$r;
        }
      }else{
        my $rev=substr($unit[0],8,1);
        if ($rev eq "f"){
           $beshit{$unit[3]}=[1,0];
        }else{
           $beshit{$unit[3]}=[0,1];
        }
      }
      if ($_=~/$ARGV[0]/){
       #my @unit=split(" ",$_);
       my $toend=$unit[10]-$unit[11];
       if (exists $ref{$unit[3]}){
           $ref{$unit[3]}.="\t$unit[5]\t$unit[10]\t$toend\t$unit[0]";
           $ctg{$unit[3]}.="\t$unit[4]";
       }else{ 
           $ref{$unit[3]}="$unit[5]\t$unit[10]\t$toend\t$unit[0]"; 
           $ctg{$unit[3]}="$unit[4]";
       }
       #print "$_\n";
      }elsif($_=~/$ARGV[1]/){
       #my @unit=split(" ",$_);
       my $toend=$unit[10]-$unit[11];
       if (exists $qry{$unit[3]}){  
           $qry{$unit[3]}.="\t$unit[5]\t$unit[10]\t$toend\t$unit[0]";
       }else{
           $qry{$unit[3]}="$unit[5]\t$unit[10]\t$toend\t$unit[0]";
       } 
       #print "$_\n";
      }
     }
}
close IN;

#open OUT, ">../output/$ARGV[0]VS$ARGV[1].beshit";
open OUT, ">>../output/scaffold.beshit" or die "$!";
foreach (keys %ref){
     if (exists $qry{$_}){
        print  "$ref{$_}\t$ctg{$_}\t$qry{$_}\t$bes{$_}\tFhit\t$beshit{$_}->[0]\tRhit\t$beshit{$_}->[1]\n";
     }
}
close OUT;
