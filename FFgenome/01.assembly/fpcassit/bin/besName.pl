#!/usr/bin/perl

open IN, "brachyantha.bes" or die "$!";
open OUT, ">brachyanthaBES.fa" or die "$!";
while (<IN>){
    chomp $_;
    if ($_=~/>OB__B(\w+)\.(\w+)/){
        my $head=$1.$2;
        print OUT ">$head\n"; 
    }else{
        print OUT "$_\n";
    }
}
close OUT;
close IN;
