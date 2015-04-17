#!/usr/bin/perl
use Getopt::Long;

my %opt;
GetOptions(\%opt,"gff:s","refFLAT:s","help");

if ($opt{help}){
   print "Usage: perl $0 -g all.gff -r refFlat.txt\n";
}

my $chr;
my $strand;
my $function;
my $gene; ### name of locus
$/="###\n";
open IN, "$opt{gff}" or die "$!";
open OUT, ">$opt{refFLAT}" or die "$!";
while (<IN>){
    my @unit=split("\n",$_);
    my $name;### name of trancript
    my $tss; ### transcription start site
    my $tes; ### transcription end site
    my $css; ### coding start site
    my $ces; ### coding end site
    my @exons; ### exon start array
    my @exone; ### exon end array
    my $exonnumber;
    foreach(@unit){
       if ($_=~/^##/){
          next;
       }else{
          my @word=split("\t",$_);
          if ($word[2] eq "gene"){
             $chr=$word[0];
             $chr=~tr/CHR/chr/;
             $strand=$word[6];
             if ($word[8]=~/Name=(.*);Alias=(.*)/){
                 $gene=$2;
                 $function=$1;
                 $function=~s/%\d+/ /g;
             }elsif($word[8]=~/ID=(.*);Name=(.*)/){
                 $gene=$1;
                 $function=$2;
                 $function=~s/%\d+/ /g;
             } 
          }elsif($word[2] eq "mRNA"){
             if ($word[8]=~/Alias=(.*)/){
                 $name=$1;
             }elsif($word[8]=~/ID=(.*);/){
                 $name=$1;
             }
             $tss =$word[3];
             $tes =$word[4];
          }elsif($word[2] eq "CDS"){
             $exonnumber++;
             push (@exons,$word[3]);
             push (@exone,$word[4]);
          }
       }
    }
    $css=$exons[0];
    $ces=$exone[$#exone];
    my $exonsline=join(",",@exons);
    my $exoneline=join(",",@exone);   
    print OUT "$gene\t$name:$function\t$chr\t$strand\t$tss\t$tes\t$css\t$ces\t$exonnumber\t$exonsline\t$exoneline\n";
}
close IN;
close OUT;
