#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","iprscan:s","project:s","help");


my $help=<<USAGE;
Get iprscan result for gene in id list.
perl $0 --table /home/jfchen/FFproject/FFgenome/03.evolution/tandem/bin/IRGSP.20.tandem.id.txt --iprscan ../input/iprscan/IRGSP.iprscan --project IRGSp
gene id have no suffix
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
    my $gene=$unit[$col];
    $hash{$gene}=1;
}
close IN;
return \%hash;
}


