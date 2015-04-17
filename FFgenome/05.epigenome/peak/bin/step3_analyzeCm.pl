#!/usr/bin/perl
=header
The script is design to analyze the raw result from sumCm.pl(BSmethyl.minus.status,BSmethyl.plus.status).
The result will be tables can be used to draw figure in Excell or other perl scripts.
=cut
use Getopt::Long;
our %opt;
GetOptions(\%opt,"ref:s","plus:s","minus:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -plus BSmethyl.plus.status -minus BSmethyl.minus.status\n";
   exit();
}

our $chrseq=fasta($opt{ref}); ### store chr sequence in hash reference $chrseq

##########################################################
my $hashplus=classCrefplus($opt{ref}); ### report C class in reference genome plus strand
my ($CG,$CHG,$CHH,$C,$mCG,$mCHG,$mCHH,$mC,$CGfreq,$CHGfreq,$CHHfreq)=methylate($opt{plus},$hashplus);
print "Plus/unMeth\t$CG\t$CHG\t$CHH\t$C\n";
print "Plus/Meth\t$mCG\t$mCHG\t$mCHH\t$mC\n";
my $freq1=join("\t",@$CGfreq);
my $freq2=join("\t",@$CHGfreq);
my $freq3=join("\t",@$CHHfreq);
print "CGfreq\t$freq1\nCHGfreq\t$freq2\nCHHfreq\t$freq3\n";
undef $hashplus;

my $hashminus=classCrefminus($opt{ref}); ### report C class in reference genome minus strand
my ($CG,$CHG,$CHH,$C,$mCG,$mCHG,$mCHH,$mC,$CGfreq,$CHGfreq,$CHHfreq)=methylate($opt{minus},$hashminus);
print "Minus/unMeth\t$CG\t$CHG\t$CHH\t$C\n";
print "Minus/Meth\t$mCG\t$mCHG\t$mCHH\t$mC\n";
my $freq1=join("\t",@$CGfreq);
my $freq2=join("\t",@$CHGfreq);
my $freq3=join("\t",@$CHHfreq);
print "CGfreq\t$freq1\nCHGfreq\t$freq2\nCHHfreq\t$freq3\n";
undef $hashminus;
#############################################################



### too slow to use. give a chr and pos to the sub, it will return you a C class, CG, CHG, CHH, C(in the end of sequence of followed by N).
sub classCpos
{
   my ($chr,$pos)=@_;
   my $seq=$chrseq->{$chr};
   my $length=length $seq;
   my $base  =substr($seq,$pos,1); 
   if ($base eq "C"){
      if ($pos == $length-1){  ### stop at last three base of sequence
          return "C";
      }elsif($pos == $length-2){
          my $first =substr($seq,$pos+1,1);
          if ($first eq "G"){
             return "CG";
          }else{
             return "C";
          } 
      }else{
         my $first =substr($seq,$pos+1,1);
         my $second=substr($seq,$pos+2,1); 
         if ($first eq "G"){ 
            return "CG";
         }elsif($second eq "G"){ 
            return "CHG";
         }elsif($first eq "N" or $second eq "N"){
            return "C";
         }else{
            return "CHH";
         }
      }
   }elsif($base eq "G"){ 
      if ($pos == 0){ ### next if in the first two base of sequence
          return "C";   
      }elsif($pos == 1){
          my $first =substr($seq,$pos-1,1);
          if ($first eq "C"){
             return "CG";
          }else{
             return "C";
          }
      }else{
          my $first =substr($seq,$pos-1,1); 
          my $second=substr($seq,$pos-2,1); 
          if ($first eq "C"){ 
             return "CG";
          }elsif($second eq "C"){ 
             return "CHG";
          }elsif($first eq "N" or $second eq "N"){ 
             return "C";
          }else{
             return "CHH";
          } 
      }
   }
}


### classify C in reference genome into CG,CHG and CHH. In minus strand,CG,CHG,HHG. 
sub classCrefplus
{
my ($file)=@_;
my %hash;
$/=">";
my ($plusCG,$plusCHG,$plusCHH);
open IN, "$file" or die "$!";
while (<IN>){
   chomp $_;
   next if (length $_ < 2);
   my @unit=split("\n",$_);
   my $head=shift @unit;
   my $seq =join("",@unit);
   $seq=~s/\t//g;
   $seq=~s/\n//g;
   $seq=~s/\s//g;
   $seq=~s/>//g;
   my $length=length $seq;
   while ($seq=~/C/ig){
      my $pos=pos($seq)-1;
      #print "$head\t$pos\n";
      next if ($pos >= $length-2);  ### stop at last three base of sequence
      my $base  =substr($seq,$pos,1);
      my $first =substr($seq,$pos+1,1);
      my $second=substr($seq,$pos+2,1);
      my $index=$head."_"."$pos";
      if ($first=~/G/i){
         $plusCG++;
         $hash{$index}="CG";
      }elsif($second=~/G/i){
         $plusCHG++;
         $hash{$index}="CHG";
      }else{
         $plusCHH++;
         $hash{$index}="CHH";
      }
   }
}
close IN;
$/="\n";
### summarize C classification in reference genome

print "Reference plus\n";
print "Class\tCG\tCHG\tCHH\n";
print "Plus\t$plusCG\t$plusCHG\t$plusCHH\n";

return \%hash;
}
##################
sub classCrefminus
{
my ($file)=@_;
my %hash;
$/=">";
my ($minusCG,$minusCHG,$minusCHH);
open IN, "$file" or die "$!";
while (<IN>){
   chomp $_;
   next if (length $_ < 2);
   my @unit=split("\n",$_);
   my $head=shift @unit;
   my $seq =join("",@unit);
   $seq=~s/\t//g;
   $seq=~s/\n//g;
   $seq=~s/\s//g;
   $seq=~s/>//g;
   my $length=length $seq;
   while ($seq=~/G/ig){
      my $pos=pos($seq)-1;
      next if ($pos <= 1); ### next if in the first two base of sequence
      my $base  =substr($seq,$pos,1);
      my $first =substr($seq,$pos-1,1);
      my $second=substr($seq,$pos-2,1);
      my $index=$head."_"."$pos";
      if ($first=~/C/i){
         $minusCG++;
         $hash{$index}="CG";
      }elsif($second=~/C/i){
         $minusCHG++;
         $hash{$index}="CHG";
      }else{
         $minusCHH++;
         $hash{$index}="CHH";
      }
   }
}
close IN;
$/="\n";

##########
print "Reference minus\n";
print "Class\tCG\tCHG\tCHH\n";
print "Minus\t$minusCG\t$minusCHG\t$minusCHH\n";

return \%hash;
}



### determine methylation status for each base
sub methylate
{
my ($file,$hash)=@_;
my (@CGfreq,@CHGfreq,@CHHfreq);
my ($CG,$CHG,$CHH,$C);
my ($mCG,$mCHG,$mCHH,$mC);
open IN, "$file" or die "$!";
open OUT1, ">$file.meth" or die "$!";
open OUT2, ">$file.unmeth" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_); 
    my $chr=$unit[0];
    my $pos=$unit[1];
    my $index=$chr."_".$pos;
    my $mClevel=$unit[4];
    if ($unit[5] eq "Me"){
        #my $Cclass=classCpos($chr,$pos);
        my $Cclass=$hash->{$index};
        if ($Cclass eq "C"){
          $mC++;
        }elsif($Cclass eq "CG"){
          $mCG++;
          if ($mClevel == 1){
             $CGfreq[9]++;
          }else{
             my $temp=int ($mClevel*10);
             $CGfreq[$temp]++;
          }
        }elsif($Cclass eq "CHG"){
          $mCHG++;
          if ($mClevel == 1){
             $CHGfreq[9]++;
          }else{
             my $temp=int ($mClevel*10);
             $CHGfreq[$temp]++;
          }
        }elsif($Cclass eq "CHH"){
          $mCHH++;
          if ($mClevel == 1){
             $CHHfreq[9]++;
          }else{
             my $temp=int ($mClevel*10);
             $CHHfreq[$temp]++;
          }
        }
        print OUT1 "$_\t$Cclass\n";
    }else{
        #my $Cclass=classCpos($chr,$pos);
        my $Cclass=$hash->{$index};
        if ($Cclass eq "C"){
          $C++;
        }elsif($Cclass eq "CG"){
          $CG++;
        }elsif($Cclass eq "CHG"){
          $CHG++;
        }elsif($Cclass eq "CHH"){
          $CHH++;
        }
        print OUT2 "$_\t$Cclass\n";
    }
}
close IN;
close OUT1;
close OUT2;
return ($CG,$CHG,$CHH,$C,$mCG,$mCHG,$mCHH,$mC,\@CGfreq,\@CHGfreq,\@CHHfreq);
}

### call methylation C based on binomial distribution, false methylation rate(non-conversion and T->C sequenceingerrors) and read depth.
sub callCm
{
my ($depth,$Cm,$cutoff)=@_;
if ($depth <= 1 or $Cm == 0){
   return 0;
}else{
  for(my $i=0;$i<=$depth;$i++){
    my $prob=pbinom($i,$depth,0.0023);
    my $k;
    if ($prob <= $cutoff){
       $k=$i-1 if $i >= 1;
       if ($Cm > $k ){
          return 1;
       }else{
          return 0;
       }
    }
  }
}
}

### return probility of binomial distribution 
sub pbinom
{
my ($x,$n,$p)=@_;
my $c;
if ($x < 1 or $x==$n){
  $c=1;
}else{
  $c=fac($n)/(fac($x)*fac($n-$x));
}
my $k1=$p**$x;
my $k2=(1-$p)**($n-$x);
my $p=$c*$k1*$k2;
return $p;
}


### return factorial of a number n: n!=n(n-1)...1
sub fac
{
my ($n)=@_;
if ($n < 2){
   return $n;
}else{
   return $n*fac($n-1);
}
}

### get system time

sub get_time {
   ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday)
     = (localtime(time));
 
   $sec  = ($sec < 10)? "0$sec":$sec;
   $min  = ($min < 10)? "0$min":$min;
   $hour = ($hour < 10)? "0$hour":$hour;
   $day  = ($day < 10)? "0$day":$day;
   $mon  = ($mon < 9)? "0".($mon+1):($mon+1);
   $year += 1900;
   
   print "$year-$mon-$day $hour:$min:$sec.00\n";
   
}


### fasta file reader
sub fasta
{
my ($file)=@_;
$/=">";
my %fastaseq;
open IN, "$file" or die "$!";
while (<IN>){
   chomp $_;
   next if (length $_ < 2);
   my @unit=split("\n",$_);
   my $head=shift @unit;
   my $seq =join("",@unit);
   $seq=~s/\t//g;
   $seq=~s/\n//g;
   $seq=~s/\s//g;
   $seq=~s/>//g;
   $fastaseq{$head}=$seq;
}
close IN;
$/="\n";
return \%fastaseq;
}
