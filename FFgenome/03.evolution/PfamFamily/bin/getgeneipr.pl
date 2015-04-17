#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","iprscan:s","project:s","help");


my $help=<<USAGE;
Get iprscan result for gene in id list.
perl $0 --table --iprscan
gene id have suffix, need to fix out. like OBR_GLEAN_10025455_OBRACH
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgene=readtable($opt{table},"0");
open OUT, ">$opt{project}.iprscan" or die "$!";
open IN, "$opt{iprscan}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $head=shift @unit;
    my $tail=join("\t",@unit);
    if (exists $refgene->{$head}){
       print OUT "$head\t$tail\n";
    }
}
close IN;
close OUT;

###########
sub readtable
{
my ($file,$col)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $gene;
    if ($unit[$col]=~/(.*)\_(\w+)/){
       $gene=$1;
    }
    $hash{$gene}=1;
}
close IN;
return \%hash;
}


