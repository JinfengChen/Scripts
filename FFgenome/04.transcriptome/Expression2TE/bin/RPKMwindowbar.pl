#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for window view RPKM gene expression level on genome.
perl $0 -chrlen IRGSP.chrlen -bar rice_RPKM_bar -gff RAP3.gff3.nr.gff.chr > log 2> log2 &
-chrlen: chr length
-bar: directory for output bar file
-gff: directory for input gff file
-rpkm: rpkm table
=cut

use Getopt::Long;
use warnings;
use strict;

our %opt;
GetOptions(\%opt,"window:s","step:s","chrlen:s","bar:s","gff:s","rpkm:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -chrlen rice.chrlen -bar rice_RPKM_bar -gff RAP3.gff3.nr.gff.chr -rpkm RAP3.rpkm\n";
   exit();
}
$opt{window} ||= 500000;
$opt{step} ||= 50000;

`mkdir $opt{bar}` unless (-e $opt{bar}); 

my $refchr=chrlen($opt{chrlen});
my $refrpkm=rpkmexpr($opt{rpkm});

my @file=glob("$opt{gff}/chr*");
chomp @file;
foreach(@file){
    my $file=$_;
    if ($_=~/(chr\d+)/){
       my $chr=$1;
       my $type="RPKM";
       print "$file\n";
       GFFdensity($file,$chr,$type,$refchr->{$chr},$opt{window},$opt{step},$refrpkm);
    }
}



#system ("/home/jfchen/software/barloader/barloader $opt{bar}");

##############
sub GFFdensity
{
my ($file,$chr,$type,$chrlen,$win,$step,$refrpkm)=@_;
my %hash;
my %meth;
open IN, "$file" or die "$!";
while (<IN>){
     chomp $_;
     #print "$_\n";
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     next if ($_ =~/^#/);
     next if ($_ =~/CDS/);
     my $line=$_;
     #print "$line\n";
     sliding($line,$win,$step,\%hash,\%meth,$refrpkm);
}
close IN;

open OUT, ">$opt{bar}/$chr.$type.bar.txt" or die "$!";
     foreach(sort {$a <=> $b} keys %meth){
          my $mean=mean($meth{$_});
          print OUT "$chr\t$_\t$mean\n";
          #my $density=100*$meth{$_}/$win;
          #print OUT "$chr\t$_\t$density\n";
     }
close OUT;
}

#################
sub sliding
{
my ($line,$win,$step,$hash,$meth,$refrpkm)=@_;
my @unit=split("\t",$line);
my $id=$1 if ($unit[8] =~/ID=(.*)?;/);
my $rpkm = $refrpkm->{$id}->[1] ? $refrpkm->{$id}->[1] : "NA";
my $bin=int ($unit[3]/$step) + 1;
for(my $i=$bin;$i>0;$i--){
   my $start=($i-1)*$step;
   my $end=$start+$win;
   #print "$i\t$unit[1]\t$start\t$end\n";
   unless ($unit[3] >= $end or $unit[4] <= $start){
        push (@{$meth->{$start}},$rpkm);
   }else{
        return;
   }
=pod
   if ($unit[3] >= $start and $unit[4] <= $end) {
        $meth->{$start}+=$unit[4]-$unit[3]+1;
   }elsif($unit[3] < $start and $unit[4] > $start){
        $meth->{$start}+=$unit[4]-$start+1;
   }elsif($unit[3] < $end and $unit[4] > $end){
        $meth->{$start}+=$end-$unit[3]+1;
   }elsif($unit[3] >= $end or $unit[4] <= $start){
        return;
   }
=cut
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


####################
sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split(" ",$_);
    if ($unit[1]=~/\-/ or $unit[1] == 0){
       #$hash{$unit[0]}=[-2,0.01];
       $hash{$unit[0]}=["NA","NA"];
    }else{
       my $temp=log10($unit[1]);
       $hash{$unit[0]}=[$temp,$unit[1]];
    }
}
close IN;
return \%hash;
}


sub mean
{
my ($num)=@_;
my $loop=0;
my $total;
foreach  (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
return 0 if ($loop == 0);
my $mean=$total/$loop;
return $mean;
}

sub log10 {
    my ($n) = shift;
    return log($n)/log(10);
}


