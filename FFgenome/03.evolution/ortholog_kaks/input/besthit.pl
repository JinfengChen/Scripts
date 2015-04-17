#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"table:s","help");


my $help=<<USAGE;
perl $0 --table
LOC_Os01g01010.1        LOC_Os01g68010.1        2e-10
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my %hash;
open IN, "$opt{table}" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   next unless ($unit[0]=~/Os/ and $unit[1]=~/Os/);
   if (exists $hash{$unit[0]}){
      $hash{$unit[0]}= $unit[2] < $hash{$unit[0]}->[1] ? [$unit[1],$unit[2]] : $hash{$unit[0]};
   }else{
      $hash{$unit[0]}= [$unit[1],$unit[2]];
   }
}
close IN;


foreach my $p (keys %hash){
   my $chr="chr".$1 if ($p=~/LOC_Os(\d+)g/);
   open OUT, ">>$chr.besthit" or die "$!";
      print OUT "$p\t$hash{$p}->[0]\t$hash{$p}->[1]\n";
   close OUT;
}

