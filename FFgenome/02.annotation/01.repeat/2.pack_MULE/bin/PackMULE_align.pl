#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","out:s","help");


my $help=<<USAGE;
perl $0 --out --fasta

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refseq=getfastaseq($opt{fasta});
& findmule($opt{out},$refseq);


sub findmule
{
my ($file)=@_;
my ($strand,$element,$start,$end);
my %type;
my %size;
my $mule;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split(" ",$_);
   #print "$unit[8]\t$unit[9]\n";
   $size{$unit[9]}->[0]++;
   $size{$unit[9]}->[1]+=abs ($unit[6]-$unit[5]);
   if (defined $strand and $unit[8] ne $strand and $unit[9] eq $element){
      ## find TSD
      $end    =$unit[6];
      my $seq =$refseq->{$unit[4]};
      my $leftP =$start-11; 
      my $rightP=$end;
      my $left  =substr($seq,$leftP,11);
      my $right =substr($seq,$rightP,11);
      my ($tsd,$align)   =findTSD($left,$right);
      if ($tsd){
         $type{1}++;
         $mule+=abs ($end-$start+1);
         print "withTSD:$start\t$end\t$element\n";
      }else{
         $type{2}++;
         print "NoTSD:$start\t$end\t$element\n";
         #print "$align";
      }
      undef $strand;
   }else{
      $strand =$unit[8];
      $element=$unit[9];
      $start  =$unit[5];
      $end    =$unit[6];
   }
}
close IN;
print STDERR "Total Length of element with TSD:$mule\n";
my ($sum1,$sum2);
foreach my $m (sort keys %size){
   print STDERR "$m\t$size{$m}->[0]\t$size{$m}->[1]\n";
   $sum1+=$size{$m}->[0];
   $sum2+=$size{$m}->[1];
}
print STDERR "Total\t$sum1\t$sum2\n";
print "WithTSD:$type{1}\nNoTSD:$type{2}\n";
}


####
sub findTSD
{
my ($left,$right)=@_;
my ($alignout,$alignlen,$identity,$seq1end,$seq2end)=align($left,$right);
#print "$alignlen\t$identity\n";
if (($alignlen == 7 or $alignlen == 8) and $identity >= $alignlen-1){
    return (1,$alignout);
}elsif($alignlen > 8 and $identity >= $alignlen-2){
    return (1,$alignout);
}else{
    return (0,$alignout);
}
}


sub getfastaseq
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
    $hash{$head}=$seq;
}
close IN;
$/="\n";
return \%hash;
}

sub writefile
{
my ($file,$line);
open OUT, ">$file" or die "$!";
     print OUT "$line";
close OUT;
}


## do local alignment using Smith-Waterman Algorithm
sub align
{
my ($seq1,$seq2)=@_;
$seq1=~tr/atgc/ATGC/;
$seq2=~tr/atgc/ATGC/;
my $match     =1;
my $mismatch  =-1;
my $gap       =-1;


## iniialization: 
my @matrix;
$matrix[0][0]{score}   =0;
$matrix[0][0]{pointer} ="none";
for (my $j =1 ;$j <= length($seq1) ;$j++){
    $matrix[0][$j]{score}  =0;
    $matrix[0][$j]{pointer}="none"; 
}
for (my $i =1 ;$i <= length($seq2) ;$i++){
    $matrix[$i][0]{score}  =0;
    $matrix[$i][0]{pointer}="none"; 
}

## fill score matrix
my $max_i  =0;
my $max_j  =0;
my $max_score =0;

for (my $i=1;$i <= length($seq2); $i++){
    for (my $j=1;$j <= length($seq1); $j++){
        my ($diagonal_score,$left_score,$up_score);
        
        ## calculate match score
        my $letter1 =substr($seq1,$j-1,1);
        my $letter2 =substr($seq2,$i-1,1);
        if ($letter1 eq $letter2){
           $diagonal_score =$matrix[$i-1][$j-1]{score}+$match;
        }else{
           $diagonal_score =$matrix[$i-1][$j-1]{score}+$mismatch;
        }

        ## calulate gap scores

        $up_score  =$matrix[$i-1][$j]{score}+$gap;
        $left_score=$matrix[$i][$j-1]{score}+$gap;
        if ($diagonal_score <= 0 and $up_score <=0 and $left_score <= 0){
              $matrix[$i][$j]{score}   =0;
              $matrix[$i][$j]{pointer} ="none";
              next; ## no need to choose best score in the next step, go to next iteration.
        }

        ## choose best score
        
        if ($diagonal_score >= $up_score){
             if ($diagonal_score >=$left_score){
                 $matrix[$i][$j]{score}  = $diagonal_score;
                 $matrix[$i][$j]{pointer}= "diagonal";
             }else{
                 $matrix[$i][$j]{score}  = $left_score;
                 $matrix[$i][$j]{pointer}= "left";
             }
        }else{
             if ($up_score >= $left_score){
                 $matrix[$i][$j]{score}  = $up_score;
                 $matrix[$i][$j]{pointer}= "up";
             }else{
                 $matrix[$i][$j]{score}  = $left_score;
                 $matrix[$i][$j]{pointer}= "left";
             }
        }
        
        ## set maximum score

        if ($matrix[$i][$j]{score} > $max_score){
             $max_i    = $i;
             $max_j    = $j;
             $max_score= $matrix[$i][$j]{score};
        }
    }
}

# trace-back

my $align1="";
my $align2="";
my $j =$max_j;
my $i =$max_i;
my $totalscore=0;
while (1){
      last if $matrix[$i][$j]{pointer} eq "none";
      $totalscore+=$matrix[$i][$j]{score};
      if ($matrix[$i][$j]{pointer} eq "diagonal"){
          $align1 .=substr($seq1, $j-1,1);
          $align2 .=substr($seq2, $i-1,1);
          $i--;
          $j--;
      }elsif($matrix[$i][$j]{pointer} eq "left"){
          $align1 .=substr($seq1, $j-1,1);
          $align2 .="-";
          $j--;
      }elsif($matrix[$i][$j]{pointer} eq "up"){
          $align1 .="-";
          $align2 .=substr($seq2, $i-1,1);
          $i--;
      }
}

$align1 = reverse $align1;
$align2 = reverse $align2;
my $alignout="$align1\n$align2\n";
#print "max_score\tmax_i\tmax_j\n";
#print "$max_score\t$max_i\t$max_j\n";
#print "Score: $totalscore\n";
#print "$align1\n"; 
#print "$align2\n";
my $seq2start=$i+1;
my $seq1start=$j+1;
my $seq2end  =$max_i;
my $seq1end  =$max_j;
#print "offset query:  $offset1\t$max_j\n";
#print "offset target: $offset2\t$max_i\n";
my $alignlen=length $align1;
my $identity;

for (my $i=0;$i < $alignlen; $i++){
     my $word1=substr($align1,$i,1);
	 my $word2=substr($align2,$i,1);
	 if ($word1 eq $word2) {
	      $identity++;
	 }
}
return ($alignout,$alignlen,$identity,$seq1end,$seq2end);
}

