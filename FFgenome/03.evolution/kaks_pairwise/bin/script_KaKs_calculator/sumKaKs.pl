#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Summary the result of KaKs_calculator and calculate_cds_aa_identity.pl.
Output informations:
Gene	Ka	Ks	Ka/Ks	Mys	CDS identity	Protein identity

Run: perl sumKaKs.pl 43890.fas.axt.kaks 43890.fas.axt.identity > sumKaKs.txt
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 
my $head;
if ($ARGV[0]=~/.*\/(.*)\.fas.axt.kaks/){
   $head=$1;
}
open IN, "$ARGV[0]" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   next if ($_ eq "");
   my @unit=split("\t");
   if ($unit[1] eq "NG"){
       my @word=split("&",$unit[0]);
       my $gene=$head;
       my $ka=$unit[2];
       my $ks=$unit[3];
       my $kaks=$unit[4];
       my $time=$unit[3]*1000/13; 
       #print "$gene\t$ka\t$ks\t$kaks\t$time";
       print "$word[0]\t$word[1]\t$ka\t$ks\t$kaks\t$time";
   }
}
close IN;

unless (-e $ARGV[1]){
   print "\n";
   exit();
}
open IN, "$ARGV[1]" or die "$!";
<IN>;
my @unit=split("\t",<IN>);
chomp @unit;
my @word=split("&",$unit[0]);
my $gene=$word[0];
my $cdsi=$unit[1];
my $proi=$unit[2];
print "\t$cdsi\t$proi\n";
close IN;
