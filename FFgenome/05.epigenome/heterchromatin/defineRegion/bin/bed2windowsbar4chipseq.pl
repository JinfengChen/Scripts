#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for sliding window view chip-seq density on genome.
perl ../bin/bed2windowsbar4chipseq.pl -chrlen IRGSP.chrlen -bar rice_chipseq_bar -type H3K4 -bed OS.H3K4.bed.chr > log 2> log2 &
-chrlen: chr length
-bar: directory for output bar file
-bed: directory for input bed file
-type: H3K4
=cut

use Getopt::Long;
use warnings;
use strict;

our %opt;
GetOptions(\%opt,"window:s","step:s","chrlen:s","bar:s","type:s","bed:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -chrlen rice.chrlen -bar rice_chipseq_bar -type H3K4 -bed OS.H3K4.bed.chr\n";
   exit();
}
$opt{window} ||= 500000;
$opt{step} ||= 50000;

`mkdir $opt{bar}` unless (-e $opt{bar}); 

my $refchr=chrlen($opt{chrlen});

my @file=glob("$opt{bed}/chr*");
chomp @file;
foreach(@file){
    my $file=$_;
    if ($_=~/(chr\d+)/){
       my $chr=$1;
       my $type=$opt{type};
       print "$file\n";
       BEDdensity($file,$chr,$type,$refchr->{$chr},$opt{window},$opt{step});
    }
}



#system ("/home/jfchen/software/barloader/barloader $opt{bar}");

##############
sub BEDdensity
{
my ($file,$chr,$type,$chrlen,$win,$step)=@_;


my %hash;
my %meth;
open IN, "$file" or die "$!";
while (<IN>){
     chomp $_;
     #print "$_\n";
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     next if ($_ =~/^#/);
     my $line=$_;
     sliding($line,$win,$step,\%hash,\%meth);
}
close IN;

open OUT, ">$opt{bar}/$chr.$type.bar.txt" or die "$!";
     foreach(sort {$a <=> $b} keys %meth){
          my $density=$meth{$_};
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
   if ($unit[1] >= $start and $unit[1] <= $end) {
        $meth->{$start}++;
   }else{
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

