#!/usr/bin/perl

### transfer TIGRrepeat to lib file that needed by repeatmasker
### >Name#class/subclass/family
### Usage: perl transferTIGR.pl TIGRrepeat OUTFILE > log &;

die "Usage: perl transferMULE.pl ricemule_dongying.fa riceMULE.fa > log &\n" if (@ARGV < 2);

$/=">";
my $class="DNA";
my $subclass="MuDR";
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[1]" or die "$!";
while(<IN>){
     chomp $_;
     next if (length $_ < 2); 
     my @unit=split("\n",$_);
     my $head=shift @unit;
     my $seq =join("\n",@unit);
     $seq=~s/\s+//g;
     print OUT ">$head#$class/$subclass\n$seq\n";
  
}
close OUT;
close IN;
