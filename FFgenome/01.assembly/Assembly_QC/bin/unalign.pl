#!/usr/bin/perl

### statistic for a alignment;
### 

use strict;
#use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;

my ($align,$format,$help);
GetOptions(
     "align:s" => \$align,
     "format:s" => \$format,
     "help:s" => \$help
);

die "Usage: perl $0 -a alignfile -f format > log &" if ($help);
die "Usage: perl $0 -a alignfile -f format > log &" unless (-f $align);
my $name=`basename $align`;
if ($name=~/(.*)\.fas\.maf/){
    $name=$1;
}

my @position;
open IN, "$align" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my $line=$_;
    if ($line=~/^a/){
       my $line1=<IN>;
       my @unit1=split(" ",$line1);
       my $refstart=$unit1[2];
       my $refend  =$unit1[2]+$unit1[3]-1;
       my $refseq  =$unit1[6];
       $refseq=~s/\-//g;
       my $reflen  =length $refseq;
       my $line2=<IN>;
       my @unit2=split(" ",$line2);
       my $qryseq  =$unit2[6];
       #print "$refstart\t$refend\t$reflen\n";
       push (@position,[$refstart,$refend]);
       while ($qryseq=~/(\-+)/g){
           my $len=length $1;
           my $pos=pos($qryseq);
           my $start=$pos-$len;
           my $end=$pos-1;
           print "$name\tintra$len\t$start\t$end\n";
      }

    }
}
close IN;
my $count;
my $last=$position[0][1]-1;
for (my $i=1;$i < @position; $i++){
    $count++;
    my $end=$position[$i][0]-1;
    print "$name\tinter$count\t$last\t$end\n" if ($end > $last);
    $last=$position[$i][1]+1;
}

