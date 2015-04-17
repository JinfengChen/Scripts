#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for window view read hit counts on genome.
Alignment of bed format should be sorted using msort (msort -k 1,n2 DNAmethly.bed > DNAmethly.bed.sorted).
=cut

use Getopt::Long;

my %opt;
GetOptions(\%opt,"window:s","bar:s","bed:s","sorting","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: run perl $0 -win 50 -bar DNAmethly.bar.txt -bed DNAmethly.bed.sorted\n";
   exit();
}
$opt{window} ||= 50;

if ($opt{sorting}){
   system ("msort -k 1,n2 $opt{bed} > $opt{bed}.sorted");
   system ("rm $opt{bed}");
   system ("mv $opt{bed}.sorted $opt{bed}");
}
my %hash;
my %meth;
my $chr="chr04";
open IN, "$opt{bed}" or die "$!";
open OUT, ">$opt{bar}" or die "$!";
while (<IN>){
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #print "$unit[1]\n";
    my $index=int $unit[1]/$opt{window};
    my $pos  =$index*$opt{window}+1;
    if (exists $hash{$pos} and $unit[0] eq $chr){
       $hash{$pos}++;
       if ($unit[3] >= $unit[4] and $unit[3] >= 2){
          $meth{$pos}++;
       }
    }elsif ($unit[0] eq $chr){
       $hash{$pos}=1;
       if ($unit[3] >= $unit[4] and $unit[3] >= 2){
          $meth{$pos}=1;
       }else{
          $meth{$pos}=0;
       }
    }else{
       foreach (sort {$a <=> $b} keys %hash){
          my $rate=$meth{$_}/$hash{$_};
          print OUT "$chr\t$_\t$rate\n";
       }
       undef %hash;
       undef %meth;
       $hash{$pos}=1;
       if ($unit[3] >= $unit[4] and $unit[3] >= 2){
          $meth{$pos}=1;
       }else{
          $meth{$pos}=0;
       }
       $chr=$unit[0];
    }
}
foreach (sort {$a <=> $b} keys %hash){
      my $rate=$meth{$_}/$hash{$_};
      print OUT "$chr\t$_\t$rate\n";
}
close IN;
close OUT;

#system ("/home/jfchen/software/barloader/barloader $opt{bar}");

