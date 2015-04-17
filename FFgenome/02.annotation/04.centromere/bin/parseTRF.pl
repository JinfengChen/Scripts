#!/usr/bin/perl
##parse Trf result(*.txt.html), report sequence in fasta format
my @seq;
my $tempseq;
open IN, "$ARGV[0]" or die "$!";
while (<IN>){
     chomp $_;
     if ($_=~/Consensus pattern/){
          $flag=1;
     }elsif($_=~/^$/ or $_=~/^Done/ and $flag==1){
          push (@seq,$tempseq);
          $flag=0;
          $tempseq="";
     }elsif($flag==1){
          $tempseq.=$_;
     }
}
close IN;

foreach (@seq){
     $counter++;
     my $len=length $_;
     if ($len>100){
        print ">$counter $len\n$_\n";
     }
}

