#!/usr/bin/perl
open IN, "$ARGV[0]" or die "$!";
while(<IN>){
    srand;
    my $index=rand(1);
    if ($index < 0.25){
       print $_;
    }
}
close IN;
