#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"os:s","ob:s","sb:s","help");


my $help=<<USAGE;
Merge gene from three seperate get related gene set.
perl merge.pl --os ../rice/4tree --ob ../brachyantha/4tree --sb ../sorghum/4tree

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @file=glob("$opt{os}/*.fa");
foreach my $file (@file){
   my %hash;
   my $type=$1 if ($file=~/$opt{os}\/(.*)\.fa$/);
   getfastaseq($file,\%hash);
   my $ob=$opt{ob}."/$type.fa";
   getfastaseq($ob,\%hash);
   my $sb=$opt{sb}."/$type.fa";
   getfastaseq($sb,\%hash);
   open OUT ,">$type.fa" or die "$!";
       foreach my $gene (keys %hash){
          print OUT ">$gene\n$hash{$gene}\n";
       }
   close OUT;
}


############################################3
sub getfastaseq
{
$/=">";
my ($file,$ref)=@_;
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
    $ref->{$head}=$seq;
}
$/="\n";
}
 
