#!/usr/bin/perl

##get sequence larger than 100kb
use Getopt::Long;

GetOptions(
         "infile:s"  => \$infile,
         "outfile:s" => \$outfile,
         "length:s"  => \$lengthcut
);


$/="\>";
open IN, "$infile" or die "can not open my file";
open OUT, ">$outfile" or die "can not open my outfile";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my @array=split(" ",$head);
    $head=$array[0];
    #print "$head\n";
    my $seq=join("",@unit);
    if (length $seq >= $lengthcut){
         print OUT ">$head\n$seq\n";
    }
}
close IN;
close OUT;
$/="\n";

