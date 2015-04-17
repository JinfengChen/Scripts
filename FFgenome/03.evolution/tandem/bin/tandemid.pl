#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Convert OB.tandem_repeat_gene.txt to OB.tandem.id.
Run: perl tandmeid.pl OB.tandem_repeat_gene.txt
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 
my %hash;
my $head;
if ($ARGV[0]=~/(\w+)\.tandem_repeat_gene\.txt/){
    $head=$1;
}

open IN, "$ARGV[0]" or die "$!";
open OUT, ">$head.tandem.id.txt" or die "$!";
while(<IN>){
    next if ($_ =~/^$/);
    my @unit=split(" ",$_);
    foreach(@unit){
       if (exists $hash{$_}){
           print "$hash{$_}\n";
       }else{
           $hash{$_}=1;
       }
    }
}

foreach(sort keys %hash){
    print OUT "$_\n";

}

close IN;
close OUT;
