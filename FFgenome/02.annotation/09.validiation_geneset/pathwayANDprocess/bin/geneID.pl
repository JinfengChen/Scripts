#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"ortholog:s","help");


my $help=<<USAGE;
Convert id from raw to gramenev1.4
perl $0 --ortholog gene.ortholog
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


my %trans2id;
open IN, "codetablev1.4" or die "$!";
while(<IN>){
     chomp $_;
     my @unit=split("\t",$_);
     $trans2id{$unit[1]}=$unit[3];
}
close IN;

open IN, "$opt{ortholog}" or die "$!";
open OUT, ">$opt{ortholog}.v1.4" or die "$!";
while (<IN>){
     chomp $_;
     my @unit=split("\t",$_);
     if ($unit[0]=~/Gene/){
        print OUT "$_\n";
     }else{
        my $id=$trans2id{$unit[5]};
        $unit[5]=$id;
        my $left=join("\t",@unit);
        print OUT "$left\n";
     }
}
close OUT;
close IN;


