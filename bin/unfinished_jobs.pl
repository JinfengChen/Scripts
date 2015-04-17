#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use File::Basename;

GetOptions (\%opt,"sh:s","cut:s","help");


my $help=<<USAGE;
perl $0 --sh run.sh --cut run.cut

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @file = glob("$opt{cut}/*.iprscan");

readshell($opt{sh}, $opt{cut});


##/opt/iprscan/5-44.0/interproscan.sh -iprlookup -goterms -T /bigdata/cjinfeng/00.RD/Annotation/Function_annotaion/Interpro/bin/A123_allpathlgv1 -f TSV -i /bigdata/cjinfeng/00.RD/Annotation/Function_annotaion/Interpro/bin/A123_allpathlgv1/A123_allpathlgv1.maker.proteins.fasta.cut/A123_allpathlgv1.maker.proteins.fasta.014 -o /bigdata/cjinfeng/00.RD/Annotation/Function_annotaion/Interpro/bin/A123_allpathlgv1/A123_allpathlgv1.maker.proteins.fasta.cut/A123_allpathlgv1.maker.proteins.fasta.014.iprscan;
sub readshell
{
my ($file, $dir)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my $line = $_;
    my $ofile = 'NA';
    if ($line=~/\-o (.*);/){
       $ofile = $dir."/".basename($1);
    }
     
    unless (-e $ofile){
       print "$line\n";
    }
}
close IN;
return \%hash;
}
 
