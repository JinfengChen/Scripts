#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blastm8:s","help");


my $help=<<USAGE;
perl $0 --blastm8
--blastm8: name.blast.m8
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

`perl /home/jfchen/FFproject/tools/solar/solar/solar.pl -a est2genome1 -f m8 -z $opt{blastm8} > exon.solar`;
`perl /home/jfchen/FFproject/tools/bin/bestAlign.pl exon.solar > exon.best.solar`;
findbadexon("exon.best.solar");


########
sub findbadexon
{
my ($file)=@_;
my %hash;
open OUT, ">gapexon.txt" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $exon=$unit[0];
    my $chr =$unit[5];
    next if $unit[9] < 2;
    my @qryblock=split(";",$unit[11]);
    my @refblock=split(";",$unit[12]);
    my @gap;
    my $flag=0;
    for(my $i=1;$i<@qryblock;$i++){
       my @ref=split(",",$refblock[$i]);
       my @ref0=split(",",$refblock[$i-1]);
       my $reflen=abs($ref[0]-$ref0[1]);
       my @qry=split(",",$qryblock[$i]);
       my @qry0=split(",",$qryblock[$i-1]);
       my $qrylen=$qry[0]-$qry0[1];
       my $start=$ref[0] > $ref0[1] ? $ref0[1] : $ref0[0];
       my $end  =$ref[0] > $ref0[1] ? $ref0[0] : $ref0[1];
       print "$exon\t$qrylen\t$reflen\n";
       if ($reflen > 200 and $qrylen < 20){
          push (@gap,"$chr:$start:$end:$reflen:$qrylen");
          $flag=1;
       }
    }
    my $hit=join("\t",@gap);
    print OUT "$exon\t$hit\n" if ($flag == 1);
}
close IN;
close OUT;
return \%hash;
}



sub getblat
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $exon=$unit[9];
    my @size=split(",",$unit[18]);
    my @qrystart=split(",",$unit[19]);
    my @refstart=split(",",$unit[20]);
    next if @size < 2;
    my @gap;
    my $flag=0;
    for(my $i=1;$i<@size;$i++){
       my $qrylen=$qrystart[$i]-($qrystart[$i-1]+$size[$i-1]);
       my $reflen=$refstart[$i]-($refstart[$i-1]+$size[$i-1]);
       if ($reflen > 100 and $qrylen < 20){
          $flag=1;
          my $start=$refstart[$i-1]+$size[$i-1];
          my $end  =$refstart[$i];
          my $chr  =$unit[13];
          push (@gap,"$chr:$start:$end");
       } 
    }
    my $hit=join("\t",@gap);
    print "$exon\t$hit\n" if ($flag == 1);
}
close IN;
return \%hash;
}

