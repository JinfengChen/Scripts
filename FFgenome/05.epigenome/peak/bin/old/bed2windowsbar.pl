#!/usr/bin/perl
=header
The script is designed to generate a bar.txt file for window view read hit counts on genome.
Alignment of bed format should be sorted using msort (msort -k 1,n2 DNAmethly.bed > DNAmethly.bed.sorted).
=cut

use Getopt::Long;

my %opt;
GetOptions(\%opt,"window:s","bar:s","bed:s","help");

if ($opt{help}){
   print "Usage: run perl $0 -win 50 -bar DNAmethly.bar.txt -bed DNAmethly.bed.sorted\n";
   exit();
}
$opt{window} ||= 50;

my %hash;
my $chr="chr1";
open IN, "$opt{bed}" or die "$!";
open OUT, ">$opt{bar}" or die "$!";
while (<IN>){
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $index=int $unit[1]/$opt{window};
    my $pos  =$index*$opt{window}+1;
    if (exists $hash{$pos} and $unit[0] eq $chr){
       $hash{$pos}++;
    }elsif ($unit[0] eq $chr){
       $hash{$pos}=1;
    }else{
       foreach (sort {$a <=> $b} keys %hash){
          print OUT "$chr\t$_\t$hash{$_}\n";
       }
       undef %hash;
       $hash{$pos}=1;
       $chr=$unit[0];
    }
}
close IN;
close OUT;

system ("/home/jfchen/software/barloader/barloader $opt{bar}");

