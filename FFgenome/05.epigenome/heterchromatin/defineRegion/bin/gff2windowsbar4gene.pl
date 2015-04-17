#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for sliding window view coding density on genome.
perl ../bin/gff2windowsbar4gene.pl -chrlen IRGSP.chrlen -bar rice_gene_bar -gff RAP3.gff3.nr.gff.chr > log 2> log2 &
-chrlen: chr length
-bar: directory for output bar file
-gff: directory for input gff file
=cut

use Getopt::Long;
use warnings;
use strict;

our %opt;
GetOptions(\%opt,"window:s","step:s","chrlen:s","bar:s","gff:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -chrlen rice.chrlen -bar rice_gene_bar -gff RAP3.gff3.nr.gff.chr\n";
   exit();
}
$opt{window} ||= 500000;
$opt{step} ||= 50000;

`mkdir $opt{bar}` unless (-e $opt{bar}); 

my $refchr=chrlen($opt{chrlen});

my @file=glob("$opt{gff}/chr*");
chomp @file;
foreach(@file){
    my $file=$_;
    if ($_=~/(chr\d+)/){
       my $chr=$1;
       my $type="GENE";
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
     next if ($_ =~/mRNA/);
     my $line=$_;
     print "$line\n";
     sliding($line,$win,$step,\%hash,\%meth);
}
close IN;

open OUT, ">$opt{bar}/$chr.$type.bar.txt" or die "$!";
     foreach(sort {$a <=> $b} keys %meth){
          my $density=100*$meth{$_}/$win;
          print OUT "$chr\t$_\t$density\n";
     }
close OUT;
}

#################
sub sliding
{
my ($line,$win,$step,$hash,$meth)=@_;
my @unit=split("\t",$line);
my $bin=int ($unit[3]/$step) + 1;
for(my $i=$bin;$i>0;$i--){
   my $start=($i-1)*$step;
   my $end=$start+$win;
   #print "$i\t$unit[1]\t$start\t$end\n";
   if ($unit[3] >= $start and $unit[4] <= $end) {
        $meth->{$start}+=$unit[4]-$unit[3]+1;
   }elsif($unit[3] < $start and $unit[4] > $start){
        $meth->{$start}+=$unit[4]-$start+1;
   }elsif($unit[3] < $end and $unit[4] > $end){
        $meth->{$start}+=$end-$unit[3]+1;
   }elsif($unit[3] >= $end or $unit[4] <= $start){
        return;
   }
}
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

