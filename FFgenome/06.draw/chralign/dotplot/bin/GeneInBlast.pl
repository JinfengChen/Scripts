#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blast:s","help");


my $help=<<USAGE;
Read blast table, count gene number of each species in blast.
perl $0 --blast 
LOC_Os01g01010.1        LOC_Os01g68010.1        2e-10
LOC_Os01g01010.1        LOC_Os02g48000.1        1e-25
LOC_Os01g01010.1        LOC_Os02g56570.1        2e-47

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

readtable($opt{blast});

sub readtable
{
my ($file)=@_;
my %hash;
my %ob;
my %os;
my %sb;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if ($unit[0]=~/Ob/i){
       $ob{$unit[0]}++;
    }elsif($unit[0]=~/Os/i){
       $os{$unit[0]}++;
    }elsif($unit[0]=~/Sb/i){
       $sb{$unit[0]}++;
    }
    if ($unit[1]=~/Ob/i){
       $ob{$unit[1]}++;
    }elsif($unit[1]=~/Os/i){
       $os{$unit[1]}++;
    }elsif($unit[1]=~/Sb/i){
       $sb{$unit[1]}++;
    }
}
close IN;
my $obn=keys %ob;
my $osn=keys %os;
my $sbn=keys %sb;
print "OB\t$obn\nOS\t$osn\nSB\t$sbn\n";
}

