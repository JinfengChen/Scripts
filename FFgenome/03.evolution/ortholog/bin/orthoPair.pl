#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Convert the table.OS-BR file to OS-BR.orth.
In OS-BR.orth, only 1:1 ortholog relation are included.

Run: perl orthoPair.pl table.OS-BR 
The output file will be OS-BR.orth.
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 

my $head;
if ($ARGV[0] =~/table\.(.*)/){
   $head=$1;
}else{
   print "Not a table file, example table.OS-BR\n";
   exit ();
}

my $total;
my $single;
my $orth=$head.".orth";
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$orth" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    $total++;
    my @unit=split("\t",$_);
    my @orthA=split(" ",$unit[2]);
    my @orthB=split(" ",$unit[3]);
    if (@orthA == 2 and @orthB == 2){
       $single++;
       print OUT "$orthA[0]\t$orthB[0]\n";
    }
}
close OUT;
close IN;

print "Total ortholog group: $total\n";
print "1:1 orthlog group: $single\n";
