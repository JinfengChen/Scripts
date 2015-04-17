#!/usr/bin/perl
=header
The script is design to summary methylation information from BED format of bowtie map result.
=cut
use Getopt::Long;
our %opt;
GetOptions(\%opt,"ref:s","1:s","2:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 SRR042638_1.fastq -2 SRR042638_2.fastq -project BSmethyl\n";
   exit();
}

## read all.fa sequence into hash %chrseq
print "Reading genome fasta file..\n";
$/=">";
my %chrseq;
open IN, "$opt{ref}" or die "$!";
while (<IN>){
   next if (length $_ <= 2);
   my @unit=split("\n",$_);
   my $head=shift @unit;
   my $seq =join("",@unit);
   $seq=~s/\t//g;
   $seq=~s/\r//g;
   $seq=~s/\n//g;
   $seq=~s/>//g;
   $chrseq{$head}=$seq;
}
close IN;
$/="\n";

## read solexa read into %read1 and %read2
print "Reading BS reads fastq file\n";
my $read1=fastq($opt{1});
my $read2=fastq($opt{2});


## read BED file and summary inf of C methylation
print "Reading BED file and summary the inf\n";
#summaryPlus("$opt{project}.T.bed");
summaryMinus("$opt{project}.A.bed");

sub summaryPlus
{
my ($file)=@_;
my %status;
open IN, "$file" or die "$!";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    my $read;
    my $pair;
    if ($unit[3]=~/(.*)\/(\d+)/){
        $read=$1;
        $pair=$2;
    }
    my $readseq;
    my $refseq;
    if ($pair==1){
       $readseq=$read1->{$read};
       $refseq =substr($chrseq{$unit[0]},$unit[1],45);
    }else{
       $readseq=$read2->{$read};
       $refseq =substr($chrseq{$unit[0]},$unit[1],45);
    }
    if ($unit[5] eq "+"){
       #print "$read\n$readseq\n$refseq\n";
    }else{
       my $revcom=reverse $readseq;
       $revcom=~tr/ATCGatcg/TAGCtagc/;
       $readseq=$revcom;
       #print "$read\n$revcom\n$refseq\n";
    }
    my $length=length $readseq;
    for(my $i=0;$i<$length;$i++){
       my $position=$unit[1]+$i;
       my $index=$unit[0]."_".$position;
       my $readbase=substr($readseq,$i,1);
       my $refbase =substr($refseq,$i,1);
       if ($refbase eq "C" and $readbase eq "C"){
             $status{$index}->[0]++;
             
       }elsif($refbase eq "C" and $readbase eq "T"){
             $status{$index}->[1]++;
             
       }
    }
}
close IN;

open OUT, ">$opt{project}.plus.status" or die "$!";
foreach(sort keys %status){
     #print "$_\n";
     my ($chr,$pos)=split("_",$_);
     my $depth=$status{$_}->[0]+$status{$_}->[1];
     my $Cm;
     if ($status{$_}->[0] eq ""){
       $Cm=0;
     }else{
       $Cm=$status{$_}->[0];
     }
     my $Me;
     if (callCm($depth,$Cm,0.01)){
        $Me="Me";
     }else{
        $Me="Un";
     }
     my $mClevel=$status{$_}->[0]/($status{$_}->[0]+$status{$_}->[1]);
     print OUT "$chr\t$pos\t$status{$_}->[0]\t$status{$_}->[1]\t$mClevel\t$Me\n";
}
close OUT;

}

###############

sub summaryMinus
{
my ($file)=@_;
my %status;
open IN, "$file" or die "$!";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\t",$_);
    my $read;
    my $pair;
    if ($unit[3]=~/(.*)\/(\d+)/){
        $read=$1;
        $pair=$2;
    }
    my $readseq;
    my $refseq;
    if ($pair==1){
       $readseq=$read1->{$read};
       $refseq =substr($chrseq{$unit[0]},$unit[1],45);
    }else{
       $readseq=$read2->{$read};
       $refseq =substr($chrseq{$unit[0]},$unit[1],45);
    }
    if ($unit[5] eq "+"){
       print "$read\n$readseq\n$refseq\n";
    }else{
       my $revcom=reverse $readseq;
       $revcom=~tr/ATCGatcg/TAGCtagc/;
       $readseq=$revcom;
       print "$read\n$revcom\n$refseq\n";
    }
    my $length=length $readseq;
    for(my $i=0;$i<$length;$i++){
       my $position=$unit[1]+$i;
       my $index=$unit[0]."_".$position;
       my $readbase=substr($readseq,$i,1);
       my $refbase =substr($refseq,$i,1);
       if ($refbase eq "G" and $readbase eq "G"){
            $status{$index}->[0]++;
       }elsif($refbase eq "G" and $readbase eq "A"){
            $status{$index}->[1]++;
       }
    }
}
close IN;

open OUT, ">$opt{project}.minus.status" or die "$!";
foreach(sort keys %status){
       my ($chr,$pos)=split("_",$_);
       my $depth=$status{$_}->[0]+$status{$_}->[1];
       my $Cm;
       if ($status{$_}->[0] eq ""){
           $Cm=0;
       }else{
           $Cm=$status{$_}->[0];
       } 
       my $Me;
       if (callCm($depth,$Cm,0.01)){
          $Me="Me";
       }else{
          $Me="Un";
       }
       my $mClevel=$status{$_}->[0]/($status{$_}->[0]+$status{$_}->[1]);
       print OUT "$chr\t$pos\t$status{$_}->[0]\t$status{$_}->[1]\t$mClevel\t$Me\n";
}
close OUT;

}

### call methylation C based on binomial distribution, false methylation rate(non-conversion and T->C sequenceing errors) and read depth.
sub callCm
{
my ($depth,$Cm,$cutoff)=@_;
if ($depth <= 1 or $Cm == 0){
   return 0;
}else{
  for(my $i=0;$i<=$depth;$i++){
    my $prob=pbinom($i,$depth,0.0023);
    my $k;
    if ($prob <= $cutoff){
       $k=$i-1 if $i >= 1;
       if ($Cm > $k ){
          return 1;
       }else{
          return 0;
       }
    }
  }
}
}

### return probility of binomial distribution 
sub pbinom
{
my ($x,$n,$p)=@_;
my $c;
if ($x < 1 or $x==$n){
  $c=1;
}else{
  $c=fac($n)/(fac($x)*fac($n-$x));
}
my $k1=$p**$x;
my $k2=(1-$p)**($n-$x);
my $p=$c*$k1*$k2;
return $p;
}


### return factorial of a number n: n!=n(n-1)...1
sub fac
{
my ($n)=@_;
if ($n < 2){
   return $n;
}else{
   return $n*fac($n-1);
}
}







##########sub fastq###############
sub fastq
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   if ($_=~/^@(.*)$/){
      my @unit=split(" ",$1);
      my $head=$unit[0];
      my $seq=<IN>;
      chomp $seq;
      $hash{$head}=$seq;
   }
}
close IN;
return \%hash;
}
