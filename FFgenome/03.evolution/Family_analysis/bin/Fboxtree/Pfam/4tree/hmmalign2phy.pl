#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"hmmalign:s","help");


my $help=<<USAGE;
perl $0 --hmmalign 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash;
open IN, "$opt{hmmalign}" or die "$1"; 
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split(" ",$_);
   $hash{$unit[0]}.=$unit[1];
}
close IN;

my @gene=keys %hash;
my $num=@gene;
my $len=length $hash{$gene[0]};
print "$num\t$len\n";
foreach(keys %hash){
   print "$_\t$hash{$_}\n";
}
