#!/usr/bin/perl

### statistic for a alignment;
### 

use strict;
#use warnings;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;

my ($align,$format,$help);
GetOptions(
     "align:s" => \$align,
     "format:s" => \$format,
     "help:s" => \$help
);

die "Usage: perl statAlign.pl -a alignfile -f format > log &" if ($help);
die "Usage: perl statAlign.pl -a alignfile -f format > log &" unless (-f $align);
my $name=`basename $align`;
if ($name=~/(.*)\.maf/){
    $name=$1;
}

my $stats=Bio::Align::DNAStatistics->new();
my $alignin=Bio::AlignIO->new(
                              -format => $format,
                              -file   => $align
                             );

my ($addlen,$addseqlen,$addn1,$addn2,$addgaplen1,$addgaplen2,$addcom,$adddiff,$addgaps,$addgapnum1,$addgapnum2);

while (my $aln =$alignin->next_aln){
     my $length=$aln->length();
     print "Length of alignment: $length\n";
     my ($newaln,$hash,$gaplen,$gapnum,$seqlen,$remove)=removecol($aln);
     my $comparable=$stats->pairwise_stats->number_of_comparable_bases($newaln);
     my $difference=$stats->pairwise_stats->number_of_differences($newaln);
     my $gaps=$stats->pairwise_stats->number_of_gaps($newaln);
     print "Number_of_comparable_bases: $comparable\n";
     print "Number_of_differences: $difference\n";
     print "Number_of_gaps: $gaps\n";
     $addseqlen+=$seqlen->{1};
     $addlen+=$length;
     $addn1 +=$hash->{1};
     $addn2 +=$hash->{2};
     $addgaplen1 +=$gaplen->{1};
     $addgaplen2 +=$gaplen->{2};
     $addgapnum1 +=$gapnum->{1};
     $addgapnum2 +=$gapnum->{2};
     $addcom+=$comparable;
     $adddiff+=$difference;
     $addgaps+=$gaps;
}

open OUT, ">>verify.txt" or die "$!";
print OUT "SeqName\tAlignLen\tNN in Seq1\tNN in Seq2\tGap in Seq1\tGap in Seq2\tCommon Site\tMismatch\tGap in Align\n";
print OUT "$name\t$addseqlen\t$addlen\t$addn1\t$addn2\t$addgaplen1\t$addgaplen2\t$addcom\t$adddiff\t$addgaps\t$addgapnum1\t$addgapnum2\n";
close OUT;
########################################333
sub removecol {
my ($aln)=@_;
my @position;
my %gaplen;
my %hash;
my %gapnum;
my $i;
my %seqlen;
foreach ($aln->each_seq()){
      my $seq=$_->seq();
     # my $len=length $seq;
     # print "$len\n";
      $i++;
      $gaplen{$i}=$seq=~tr/-/-/;
      $seqlen{$i}=$seq=~tr/ATCGatcg/ATCGatcg/;
      print "Gap length in Sequence $i: $gaplen{$i}\n";
      while ($seq=~/(N+)/g){
           my $len=length $1;
           $hash{$i}+=$len;
           #print "$len\n";
           my $pos=pos($seq);
           #print "$pos\n";
           my $start=$pos-$len;
           my $end=$pos-1;
           push (@position,[$start,$end]); 
      }
      
}
foreach (keys %hash){
    print "Number of N in $_ : $hash{$_}\n";
}
@position=overlap(@position);
my $remove=0;
foreach (@position){
      my $a1=$_->[0];
      my $b1=$_->[1];
      #print "$a1\t$b1\n";
      $remove+=$b1-$a1+1; 
}
print "Number of N removed in alignment : $remove\n";
my $newaln;
if (@position > 1){ 
  $newaln=$aln->remove_columns(@position);
}else{
  $newaln=$aln;
}

my $j;
foreach ($newaln->each_seq()){
      my $seq=$_->seq();
      $j++;
      while($seq=~/(-+)/g){
           $gapnum{$j}+=1;
      }
}

return ($newaln,\%hash,\%gaplen,\%gapnum,\%seqlen,$remove);
}
###############################


### merge overlap columns block by hash, 1 4, 2 7, 12 14 => 1 7, 12 14
sub overlap{
my (@array)=@_;
my %hash;
foreach (@array){
      my $a1=$_->[0];
      my $b1=$_->[1];
      for(my $i=$a1;$i<=$b1;$i++){
         $hash{$i}=1;
      }

}

my @unit=sort {$a <=> $b} keys %hash;
my @newarray;
my $start=$unit[0];
my $end=$unit[0];
for(my $i=1;$i<@unit;$i++){
     if ($unit[$i] == $end+1){
         $end=$unit[$i];
     }else{ 
         push (@newarray,[$start,$end]);
         $start=$unit[$i];
         $end=$unit[$i];   
     } 
} 
push (@newarray,[$start,$end]);
return @newarray;
}
###################################


###############################
sub show{

my ($newaln)=@_;
foreach ($newaln->each_seq()){
        my $seq=$_->seq();
        print "Show Align Seq:\n";
        print "$seq\n";
}
}
######################33
=pod
my $comparable=$stats->pairwise_stats->number_of_comparable_bases($aln);
my $difference=$stats->pairwise_stats->number_of_differences($aln);
my $gaps=$stats->pairwise_stats->number_of_gaps($aln);
print "number_of_comparable_bases: $comparable\n";
print "number_of_differences: $difference\n";
print "number_of_gaps: $gaps\n";

my $comparable=$stats->pairwise_stats->number_of_comparable_bases($newaln);
my $difference=$stats->pairwise_stats->number_of_differences($newaln);
my $gaps=$stats->pairwise_stats->number_of_gaps($newaln);
print "Number_of_comparable_bases: $comparable\n";
print "Number_of_differences: $difference\n";
print "Number_of_gaps: $gaps\n";
=cut

