#!/usr/bin/perl
use Getopt::Long;
#use warnings;
#use strict;
GetOptions (\%opt,"ortholog:s","cds1:s","cds2:s","help");


my $help=<<USAGE;
perl $0 --ortholog --cds1 --cds2
Gene1	Gene2    Ka      Ks      Ka/Ks   Myr     CDS identity    Protein identity
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

##
my $kakspipe="/home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/pipe_kaks_file.pl";
my $dir ="./kaks";
my $orth=ortholog($opt{ortholog});
my $cds1=getfastaseq($opt{cds1});
my $cds2=getfastaseq($opt{cds2});
`mkdir $dir` unless (-d $dir);
foreach my $rank (sort keys %$orth){
    my $g1=$orth->{$rank}->[0];
    my $g2=$orth->{$rank}->[1];
    my $seq1=$cds1->{$g1};
    my $seq2=$cds2->{$g2};
    my $file="$dir/$g1"."vs"."$g2".".cds.fa";
    my $line=">$g1\n$seq1\n>$g2\n$seq2\n";
    writefile($file,$line);
    `perl $kakspipe $file >> ortholog.inf`;
}

######
sub writefile
{
my ($file,$line)=@_;
open OUT, ">$file" or die "$!";
     print OUT "$line\n";
close OUT;
}

sub ortholog
{
my ($file)=@_;
my %hash;
my $count;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    $count++;
    my @unit=split("\t",$_);
    $hash{$count}=[$unit[0],$unit[1]];
}
close IN;
return \%hash;
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

