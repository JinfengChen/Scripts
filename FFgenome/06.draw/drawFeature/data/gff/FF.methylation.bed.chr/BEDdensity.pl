#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","help");

my $help=<<USAGE;
Run in where the bed chr file is.
perl $0 -fasta ffseq.fa
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 



my $win=1000000; ##1Mb
my $step=200000; ##200kb
my $reflen=getfastalen($opt{fasta});

my @file=glob("chr*");
chomp @file;
foreach(@file){
    if ($_=~/(chr\d+)/){
       my $chr=$1;
       BEDdensity($chr,$reflen->{$chr},$win,$step);
    }
}



############sub functions#####################


sub BEDdensity
{
my ($chr,$chrlen,$win,$step)=@_;
print "$chr\n";

my @bin;
my $run=int(($chrlen-$win)/$step)+1;
for (my $i=0;$i<=$run;$i++){
    my $start=$i*$step;
    my $end  =$start+$win;
    push (@bin,[$start,$end]); 
}

my %start;
open IN, "$chr" or die "$!";
while (<IN>){
     chomp $_;
     next if ($_ eq "" );
     next if ($_ =~/^chrUN/);
     my @unit=split("\t", $_);
     foreach(@bin){
         if ($unit[1] >= $_->[0] and $unit[1] <= $_->[1]){
             $start{$_->[0]}+=1;
         }
     }

}
close IN;

open OUT, ">$chr.bed.density" or die "$!";
     foreach(sort {$a <=> $b} keys %start){
          my $density=100*$start{$_}/$win;
          print OUT "$_\t$density\n";
     }
close OUT;
}



########################################################


sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=length $seq;
}
close IN;
$/="\n";
return \%hash;
}

