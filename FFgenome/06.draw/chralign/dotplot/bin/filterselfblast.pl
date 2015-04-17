#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blast:s","help");


my $help=<<USAGE;
perl $0 --blast 3wayv1.4.blast > 3wayv1.4.blast.noself
filer self blast. For example you have OS-OB-SB three species blast, we kick out these OS-OS, OB-OB, SB-SB pairs
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable($opt{blast});
##

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if (substr($unit[0],0,2) ne substr($unit[1],0,2)){
       print "$unit[0]\t$unit[1]\t$unit[2]\n";
    } 
}
close IN;
return \%hash;
}
 
