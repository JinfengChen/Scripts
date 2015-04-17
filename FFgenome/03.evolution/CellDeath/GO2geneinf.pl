#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"ipr:s","go:s","table:s","id:s","help");


my $help=<<USAGE;
Find go and iprscan annotation for list of gene from GO id.

Run: perl GO2geneinf.pl -ipr FF.ipr -go FF.GO -table FF.result -id GO:0008219
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $refipr=parseann($opt{ipr});
my $refgo =parseann($opt{go});

open IN, "$opt{table}" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    #print "$unit[1]\t$unit[8]\n";
    if ($unit[1] eq $opt{id}){
        my @gene=split(",",$unit[8]);
        foreach(@gene){
           #print "$_\n";
           if (exists $refipr->{$_}){
              print "$_\t$refipr->{$_}\n";
           }
           if (exists $refgo->{$_}){
              print "$_\t$refgo->{$_}\n";
           } 
        }
    }
}
close IN; 



##########sub functions##################
sub parseann
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my ($gene,$num,$ann)=split("\t",$_,3);
   $hash{$gene}=$ann;
}
close IN;
return \%hash;
}


