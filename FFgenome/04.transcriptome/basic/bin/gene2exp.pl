#!/usr/bin/perl

use Getopt::Long;

GetOptions(
      "soap:s"   => \$soap,
      "length:s" => \$length,
      "reads:s"  => \$reads,
      "help"     => \$help
);

if ($help){
    print "Usage: perl gene2exp.pl -s Root.soapout -l genelen -r 13452494> geneexpression\n";
    exit;
}

my %len;
open IN, "$length" or die "$!";
while (<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    $len{$unit[0]}=$unit[1];
}
close IN;


my %hash;
open IN, "$soap" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if ($unit[3] ==1 and $unit[9] ==0){
        if (exists $hash{$unit[7]}){
           $hash{$unit[7]}++;
        }else{
           $hash{$unit[7]}=1;
        }  

    }

}
close IN;
print "Gene\tLength\tReads\tRPKM\n";
foreach (sort keys %len){
   if (exists $hash{$_}){
      my $rpkm=$hash{$_}*1000000000/($len{$_}*$reads);
      print "$_\t$len{$_}\t$hash{$_}\t$rpkm\n";
   }else{
      print "$_\t$len{$_}\t0\t0\n";
   }
}

