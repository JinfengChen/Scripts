#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"inversion:s","refseq:s","qryseq:s","flanking:s","help");


my $help=<<USAGE;
perl $0 --inversion --refseq --qryseq 
Analysis flanking sequence of inversion to find invert repeat.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{flanking} ||= 10000;
my $refseq=getfastaseq($opt{refseq});
my $qryseq=getfastaseq($opt{qryseq});

my $rank;
open IN, "$opt{inversion}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   $rank++;
   my @unit=split("\t",$_);
   my $ref=invert($refseq->{"chr$unit[0]"},$unit[5],$unit[6],"ref$rank");
   my $qry=invert($qryseq->{"chr$unit[0]"},$unit[7],$unit[8],"qry$rank");
 
}
close IN;



#######################
sub invert
{
my ($seq,$temp1,$temp2,$title)=@_;
my $start=$temp1 > $temp2 ? $temp2 : $temp1;
my $end  =$temp1 > $temp2 ? $temp1 : $temp2;
my $s1=$start-$opt{flanking};
my $subseq1=substr($seq,$s1,$opt{flanking});
writefasta("$title.up.fa",$subseq1);
my $s2=$end;
my $subseq2=substr($seq,$s2,$opt{flanking});
writefasta("$title.down.fa",$subseq2);
`ssearch36 $title.up.fa $title.down.fa >> ssearch.out`;
`rm *.fa`;
}

sub ssearch
{
my ($file)=@_;
my @result;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/ or $_=~/^#/);
   push (@result,$_);
}
close IN;
return \@result;
}


sub writefasta
{
my ($file,$seq)=@_;
my $name=$1 if ($file=~/(.*)\.fa/);
open OUT, ">$file" or die "$!";
     print OUT ">$name\n$seq\n";
close OUT;
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
close IN;
$/="\n";
return \%hash;
}


