#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"id:s","help");


my $help=<<USAGE;
perl $0 --id 
Select 50 small gap genes and 50 large gap genes (> 100 bp gap).
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my %hash;
open IN, "$opt{id}" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   push (@{$hash{$unit[0]}},$unit[1]);
}
close IN;


my %large;
my %small;
foreach my $gene (keys %hash){
   next if (@{$hash{$gene}} > 1);
   next unless ($gene=~/Ob\d{2}g/);
   if ($hash{$gene}[0] > 100 and $large < 50){
      $large{$gene}=$hash{$gene}[0];
      $large++;
   }
   if ($hash{$gene}[0] < 100 and $small < 50){
      $small{$gene}=$hash{$gene}[0];
      $small++;
   }
}
print "";
foreach (sort keys %large){
   print "1\t$_\t$large{$_}\n";
}
foreach (sort keys %small){
   print "2\t$_\t$small{$_}\n";
}


