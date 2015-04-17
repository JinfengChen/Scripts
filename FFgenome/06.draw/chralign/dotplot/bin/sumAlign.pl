#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"align:s","bed:s","help");


my $help=<<USAGE;
Summary mcscan align file.
Report 1, length of genomic interval between two colinear genes (from this gene to last gene).
chr     gene1                   start   length  gene2                   start   length
Os01    LOC_Os01g01030.1        11648   1831    OBR_GLEAN_10003597      47693   25149 
Report 2, length of genomic interval in 1 Mb region of rice.
chr     bin(Mb) length of Os-Ob
Os01    33      21801 
perl $0 --align --bed

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

parsealign($opt{align});

sub parsealign
{
my ($align)=@_;
my $chr=$1 if ($align=~/chr(\d+)\.aligns/);
my $refbed;
if ($opt{bed}){
   $refbed=readbed($opt{bed});
}else{
   my $bed=$1 if ($align=~/(.*)\.aligns/);
   $refbed=readbed("$bed.bed");
}
my %hash;
my $rank;
my $flag;
open IN, "$align" or die "$!";
my $i=1;
while($i<=11){
$i++;
<IN>;
}
while(<IN>){
    chomp $_;
    if ($_=~/## Alignment (\d+)\:.*N=(\d+) Os$chr&Ob$chr plus/ and $2 > 50){
        #print "$_\n";
        $rank=$1;
        $flag=1;
    }elsif($_=~/## Alignment/){
        $flag=0;
    }elsif($flag == 1){
        my @unit=split("\t",$_); 
        push (@{$hash{$rank}},[$unit[1],$unit[2]]);
    }
}
close IN;
####
my %size;
open OUT, ">$align.table" or die "$!";
foreach my $a (keys %hash){
     my $last1;
     my $last2;
     my @pair=@{$hash{$a}};
     $last1=$refbed->{$pair[0][0]}->[1];
     $last2=$refbed->{$pair[0][1]}->[1];
     for(my $i=1;$i<@pair;$i++){
        my $interval1=$refbed->{$pair[$i][0]}->[0]-$last1 > 0 ? $refbed->{$pair[$i][0]}->[0]-$last1 : 100;
        my $interval2=$refbed->{$pair[$i][1]}->[0]-$last2 > 0 ? $refbed->{$pair[$i][1]}->[0]-$last2 : 100;
        print OUT "Os$chr\t$pair[$i][0]\t$refbed->{$pair[$i][0]}->[0]\t$interval1\t$pair[$i][1]\t$refbed->{$pair[$i][1]}->[0]\t$interval2\n";
        $last1=$refbed->{$pair[$i][0]}->[1];
        $last2=$refbed->{$pair[$i][1]}->[1];

        ##
        my $bin=1+int $refbed->{$pair[$i][0]}->[0]/1000000;
        $size{$bin}+=$interval1-$interval2;
     } 
}
close OUT;
open OUT, ">$align.size" or die "$!";
     foreach my $b (keys %size){
        print OUT "Os$chr\t$b\t$size{$b}\n";
     }
close OUT;
}


sub readbed
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[3]}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
}
 
