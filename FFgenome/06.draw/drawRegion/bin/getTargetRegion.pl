#!/usr/bin/perl

use Getopt::Long;

GetOptions (\%opt,"fasta:s","table:s","gene:s","te:s","project:s","help");


my $help=<<USAGE;
This script is designed to get sequence fregment and annotation from fasta and gff file according to position in table file.
Example table:
Query                           Target  Start   End     Strand    Length  HitLength
OB_Scaffold000110_425519_639420 chr08   4708251 5259086 +         550836  70338

Output a CompareRegion.txt file:
Example CompareRegion.txt:
OB_Scaffold000002_1830930_2143221       Os_chr09_18409023_18805233

Run: perl getTargetRegion.pl -fasta genome.fa -table HitRegionInf.txt -gene gene.gff -te te.gff -project Os

USAGE

if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit ();
} 

my $refseq=getfastaseq($opt{fasta});

open IN, "$opt{table}" or die "$!";
open OUT, ">$opt{project}_HitRegion.fa" or die "$!";
open OUT2, ">CompareRegion.txt" or die "$!";
while(<IN>){
    my @unit=split("\t",$_);
    my $len=$unit[3]-$unit[2]+1;
    my $strand=$unit[4];
    my $seq=substr($refseq->{$unit[1]},$unit[2],$len);
    if ($strand eq "-"){
       $seq=revcom($seq);  
    }
    my $head=$opt{project}."_".$unit[1]."_".$unit[2]."_".$unit[3];
    my $run1=getsubGFF($opt{gene},$head,$strand);
    my $run2=getsubGFF($opt{te},$head,$strand);
    print OUT ">$head\n$seq\n";
    print OUT2 "$unit[0]\t$head\n"; 
}
close IN;
close OUT;
close OUT2;
####################sub function############################

sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}



sub getsubGFF
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
#print "$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
       $unit[3]=$unit[3]-$start;
       $unit[4]=$unit[4]-$start;
       my $line=join("\t",@unit);
       print OUT1 "$line";
     } 

   }
   close OUT1;
   close IN1;
}else{
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
        my $tempend   =$len-($unit[3]-$start);
        my $tempstart =$len-($unit[4]-$start); 
        $unit[3]=$tempstart;
        $unit[4]=$tempend;
        if ($unit[6] eq "+"){
            $unit[6] = "-";
        }else{
            $unit[6] = "+";
        }
        my $line=join("\t",@unit);
        print OUT1 "$line";
     }
   }
   close OUT1;
   close IN1;
}
return 1;
}

sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}

