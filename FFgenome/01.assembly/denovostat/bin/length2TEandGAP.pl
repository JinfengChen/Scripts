#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"masked:s","help");


my $help=<<USAGE;
Calculate TE% and GAP% for each scaffold.
GAP are NN in sequence.
TE are masked to lower case in sequence.
Run: perl length2TEandGAP.pl -masked scaffold.masked 
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

open OUT, ">length2gapTE" or die "$!";
print OUT "ID\tLength\tGAP%\tGAP Length\tGAP Number\tTE%\tTE length\n";
my $refseq=getfastaseq($opt{masked});
foreach(sort keys %$refseq){
    my $sequence=$refseq->{$_};
    my $len=length $sequence;
    my $kb=$len/1000;
    $len=sprintf("%.2f",$kb);
    my ($gap,$gaplen,$gapnum)=estimateGAP($sequence);
    my ($te, $telen) =estimateTE($sequence);
    my $gaplenkb=$gaplen/1000;
    my $telenkb =$telen/1000;
    $gaplen=sprintf("%.2f",$gaplenkb);
    $telen =sprintf("%.2f",$telenkb);
    print OUT "$_\t$len\t$gap\t$gaplen\t$gapnum\t$te\t$telen\n";
}
close OUT;

#####################################################################3

sub estimateGAP{
my ($seq)=@_;
my $length=length $seq;
my $gapn=0;
my $gapl=0;
while ($seq=~/(N+)/g){
       $gapn++;
       $gapl+=length $1;
}

#my $n=$seq=~tr/Nn/Nn/;
#my $gapfrq=sprintf("%.2f",$n/$length);
my $gapfrq=sprintf("%.2f",$gapl/$length);
return ($gapfrq,$gapl,$gapn);
}

sub estimateTE{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
my $n=$seq=~tr/Nn/Nn/;
my $len=$length-$n;
my $te=$seq=~tr/atcg/atcg/;
if ($len > 0){
  my $tefrq=sprintf("%.2f",$te/$len);
  return ($tefrq,$te);
}else{
  my $tefrq="Na";
  return ($tefrq,0);
}
}



sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}

