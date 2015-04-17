#!/usr/bin/perl

use strict;
use FindBin qw ($Bin);
use Getopt::Long;
use SVG;

#################################Usage ##############################################
my %opt;
GetOptions(\%opt,"table:s","help:s");

my $help=<<USAGE;
perl $0 --table

USAGE

if(defined $opt{help} or keys %opt < 1){
        die  $help ;
}

my $qrychr=readtable("../input/gramenev1.4.chrlen");
my $qrycent=readtable("../input/gramenev1.4.cent");
my $refchr=readtable("../input/tigr6.chrlen");
my $refcent=readtable("../input/tigr6.cent");


my $svg=SVG->new(width=>2500,height=>800);
my $startw=100; 
my $endw  =2400;
my $starth=100;
my $endh  =700;

my $scale=getscale($refchr,"600");


#### chromosome line
my %refpos;
my %qrypos;
my @chr=sort keys %$refchr;
for (my $i=0;$i<@chr;$i++){
    ###################reference chromosome
    my $refx=$i*200+$startw;
    my $refy=$starth;
    my $refw=20;
    my $refh=$refchr->{$chr[$i]}/$scale;
    $refpos{$chr[$i]}=$refx;
    #print "$i\t$chr[$i]\t$refx\n";
    my $refline=$svg->rectangle(
              x=>$refx,y=>$refy,
              width=>$refw,height=>$refh,
              rx=>20,ry=>20,
              style=>{
                   fill=>'none',
                   stroke=>'#7A8B8B'
              }
              
    );
    my $refcx=$refx+$refw/2;
    my $refcy=$refcent->{$chr[$i]}/$scale+$refy;
    my $refrx=$refw/2;
    my $refry=$refrx;
    $svg=centromere($svg,$refcx,$refcy,$refrx,$refry,"black");
    #################query chromosome
    my $qryx=$refx+80;
    my $qryy=$refy;
    my $qryw=$refw;
    my $qryh=$qrychr->{$chr[$i]}/$scale;
    $qrypos{$chr[$i]}=$qryx;
    my $qryline=$svg->rectangle(
              x=>$qryx,y=>$qryy,
              width=>$qryw,height=>$qryh,
              rx=>30,ry=>30,
              style=>{
                   fill=>'none',
                   stroke=>'#7A8B8B'
              }
    );
    my $qrycx=$qryx+$qryw/2;
    my $qrycy=$qrycent->{$chr[$i]}/$scale+$qryy;
    my $qryrx=$qryw/2;
    my $qryry=$qryrx;
    $svg=centromere($svg,$qrycx,$qrycy,$qryrx,$qryry,"black");
}


=pod
my $refstarth=$starth+200;
my $refwidth =$endw-$startw;
my $refline  =$svg->rectangle(
              x=>$startw,y=>$refstarth,
              width=>$refwidth,height=>10,
              rx=>5.2,ry=>2.4,
              style=>{
                  stroke=>'black'
              }
);
my $refnote  =$svg->text(
              x=>60, y=>$refstarth+8,
              style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
              }
)->cdata("test");
=cut
writesvg("test.svg",$svg);

#####

sub getscale
{
my ($chr,$height)=@_;
my @len=sort {$b <=> $a} values %$chr;
#print "$len[0]\n";
my $rate=$len[0]/$height;
return $rate;
}


sub centromere
{
my ($svg,$cx,$cy,$rx,$ry,$fillcolor)=@_;
my $tag = $svg->ellipse(
        cx=>$cx, cy=>$cy,
        rx=>$rx, ry=>$ry,
        style=>{
            'fill'=>$fillcolor,
        }
    );
return $svg;
}


sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}



################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
system "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
}


