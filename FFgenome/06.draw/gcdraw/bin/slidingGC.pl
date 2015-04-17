#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;
our %opt;
GetOptions(\%opt,"help");

my $help=<<USAGE;
This program is designed to caculate gc contain of a fasta sequence in sliding windows. The results is in a file and in a dir of GCcontain that seperate every record in one file
Run: perl slidingGC.pl sample.fa
USAGE

if ($opt{help}){
   print "$help\n";

}



$/=">";
open OUT, ">$ARGV[0]\.gc" or die "can not open gc out file";
#open IN, "EC611973.txt";
open IN, "$ARGV[0]" or die "can not open my infile";
while (<IN>){

     my @unit=split("\n");
     my $temp=shift @unit;
     my @temp=split(" ",$temp);
     my $head=shift @temp;
     next if (length $head < 2);
     my $seq=join("",@unit);
     $seq=~s/\>//;
        my $length=length $seq;
        my $step=100;
        my $win=200;
        my $run=int (($length-$win)/$step)+1;
        my $gcdir="./GCcontent";
        system ("mkdir $gcdir") unless(-e $gcdir);
        open OUT1, ">./GCcontent/$head\.gc" or die "$!";
        for (my $i=0;$i<=$run;$i++){
            my $start =$i*$step;
            #print "$start\n";
            my $subseq=substr($seq,$start,$win);
           
            my $mid=$start+$win/2;
            my $gc=estimateGC($subseq);
            unless ($gc eq "Na"){
               print OUT "$mid\t$gc\n";
               print OUT1 "$mid\t$gc\n";
            }
        }
        close OUT1;
}
close IN;
close OUT;

sub estimateGC{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
my $n=$seq=~tr/Nn/Nn/;
my $len=$length-$n;
my $c=$seq=~tr/Cc/Cc/;
my $g=$seq=~tr/Gg/Gc/;
if ($len > 0){
  my $gc=($g+$c)/$len;
  return $gc;
}else{
  my $gc="Na";
  return $gc;
}
}
