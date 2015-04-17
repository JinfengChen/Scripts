#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","id:s","help");


my $help=<<USAGE;

perl $0 -table IRGSP.specific.fa.blast.table -id IRGSP.specific.id > IRGSP.specific.ann
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 


my %hit;
open IN, "$opt{table}" or die "$!";
while (<IN>){
      chomp $_;
      next if ($_ eq "");
      next if ($_ =~/Query_id/);
      my @unit=split("\t",$_);
      my $ann=pop @unit;
      next if (exists $hit{$unit[0]});
      next if ($ann =~/^Os/);
      if ($ann=~/(.*)?\[/){
          $hit{$unit[0]}=$1;
      }else{
          $hit{$unit[0]}=$ann;
      }
}
close IN;

open IN, "$opt{id}" or die "$!";
while(<IN>){
    chomp $_;
    if (exists $hit{$_}){
        print "$_\t$hit{$_}\n";
    }else{
        print "$_\tNA\n";
    }
}
close IN;
