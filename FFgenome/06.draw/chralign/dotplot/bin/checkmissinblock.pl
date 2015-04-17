#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","blockgene:s","blast:s","help");


my $help=<<USAGE;
perl $0 --bed --blockgene --blast

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $refblock=readblockgene($opt{blockgene});
my $refbed=readbed($opt{bed});
my $refblast=readblast($opt{blast});
foreach my $g (sort keys %$refbed){
   unless (exists $refblock->{$g}){
      print "$g\t$refblast->{$g}\n";
   }
}

sub readblockgene
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
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
    $hash{$unit[3]}=1;
}
close IN;
return \%hash;
}

sub readblast
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    next unless ($unit[0]=~/Os/);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

