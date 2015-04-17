#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for sliding window view methylation level on genome.
perl ../bin/bed2windowsbar4methylation_2.pl -chrlen IRGSP.chrlen -bar rice_methylation_bar -bed rice_methylation_bed
-chrlen: chr length
-bar: directory for output bar file
-bed: directory for input bed file
=cut

use Getopt::Long;
use warnings;
use strict;

our %opt;
GetOptions(\%opt,"window:s","step:s","chrlen:s","bar:s","bed:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -chrlen rice.chrlen -bar rice_methylation_bar -bed rice_methylation_bed\n";
   exit();
}
$opt{window} ||= 500000;
$opt{step} ||= 50000;

`mkdir $opt{bar}` unless (-e $opt{bar}); 

my $refchr=chrlen($opt{chrlen});

my @file=glob("$opt{bed}/*.bed");
chomp @file;
foreach(@file){
    my $file=$_;
    if ($_=~/(chr\d+)\.(\w+)\.bed/){
       my $chr=$1;
       my $type=$2;
       print "$file\n";
       BEDdensity($file,$chr,$type,$refchr->{$chr},$opt{window},$opt{step});
    }
}



#system ("/home/jfchen/software/barloader/barloader $opt{bar}");

##############
sub BEDdensity
{
my ($file,$chr,$type,$chrlen,$win,$step)=@_;

=pod
my @bin;
my $run=int(($chrlen-$win)/$step)+1;
for (my $i=0;$i<=$run;$i++){
    my $start=$i*$step;
    my $end  =$start+$win;
    push (@bin,[$start,$end]);
}
=cut

my %hash;
my %meth;
open IN, "$file" or die "$!";
while (<IN>){
     chomp $_;
     #print "$_\n";
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     my $line=$_;
=pod
     my @unit=split("\t", $_);
     foreach(@bin){
         if ($unit[1] >= $_->[0] and $unit[1] <= $_->[1]){
            $hash{$_->[0]}++;
            if ($unit[3] >= $unit[4] and $unit[3] >= 2){  ### if this loci is methylated
               $meth{$_->[0]}++;
            }
         }
     }
=cut
    sliding($line,$win,$step,\%hash,\%meth);
}
close IN;

open OUT, ">$opt{bar}/$chr.$type.bar.txt" or die "$!";
     foreach(sort {$a <=> $b} keys %hash){
          my $density=100*$meth{$_}/$hash{$_};
          print OUT "$chr\t$_\t$density\n";
     }
close OUT;
}

#################
sub sliding
{
my ($line,$win,$step,$hash,$meth)=@_;
my @unit=split("\t",$line);
my $bin=int ($unit[1]/$step) + 1;
for(my $i=$bin;$i>0;$i--){
   my $start=($i-1)*$step;
   my $end=$start+$win;
   #print "$i\t$unit[1]\t$start\t$end\n";
   if ($unit[1] >= $start and $unit[1] < $end) {
            $hash->{$start}++;
            #my @m=split(":",$unit[3]);
            #if ($m[0]+$m[1] >= 3){
            #   $hash->{$start}++;
            #   if ($m[0] >= $unit[4]){
            #       $meth->{$start}++;
            #   }
            #}
            if ($unit[3] >= $unit[4]){
            #if ($unit[3] >= $unit[4] and $unit[3] >= 2){  ### if this loci is methylated
               $meth->{$start}++;
            }
   }elsif($unit[1] >= $end or $unit[1] < $start){
            return;
   }
}
return;
}


#################
### get chr length hash
sub chrlen
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

