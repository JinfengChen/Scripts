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
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split(" ",$_);
   #print "$unit[8]\t$unit[9]\n";
   if (defined $strand and $unit[8] ne $strand and $unit[9] eq $element){
      ## find TSD
      $end    =$unit[6];
      my $seq =$refseq->{$unit[4]};
      #print "$start\t$end\n";
      #my $tsd1=findTSD1($seq,$start,$end,"11");
      #my $tsd2=findTSD1($seq,$start,$end,"10");
      #my $tsd3=findTSD1($seq,$start,$end,"9"); 
      #my $tsd4=findTSD1($seq,$start,$end,"8");
      #my $tsd5=findTSD1($seq,$start,$end,"7");
      if (findTSD1($seq,$start,$end,"11")){
         $type{1}++;
         print "withTSD:$start\t$end\t$element\n";
      }elsif(findTSD1($seq,$start,$end,"10")){
         $type{1}++;
         print "withTSD:$start\t$end\t$element\n";
      }elsif(findTSD1($seq,$start,$end,"9")){
         $type{1}++;
         print "withTSD:$start\t$end\t$element\n";
      }elsif(findTSD1($seq,$start,$end,"8")){
         $type{1}++;
         print "withTSD:$start\t$end\t$element\n";
      }elsif(findTSD1($seq,$start,$end,"7")){
         $type{1}++;
         print "withTSD:$start\t$end\t$element\n";
      }else{
         $type{2}++;
         print "NoTSD:$start\t$end\t$element\n";
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
print "WithTSD:$type{1}\nNoTSD:$type{2}\n";
}


####
sub findTSD
{
my ($left,$right)=@_;
my ($alignout,$alignlen,$identity,$seq1end,$seq2end)=align($left,$right);
print "$alignlen\t$identity\n";
if (($alignlen == 7 or $alignlen == 8) and $identity >= $alignlen-1){
    return (1,$alignout);
}elsif($alignlen > 8 and $identity >= $alignlen-2){
    return (1,$alignout);
}else{
    return (0,$alignout);
}
}
#### use emboss
sub findTSD1
{
my ($seq,$start,$end,$len)=@_;
#print "Length:$len\n";
my $leftP =$start-$len;
my $rightP=$end+1;
my $left  =substr($seq,$leftP,$len);
my $right =substr($seq,$rightP,$len);
my ($identity,$length)=align1($left,$right);
if ($len > 8){
   if ($identity >= $len-2){
      return 1;
   }else{
      return 0;
   }
}else{
   if ($identity >= $len-1){
      return 1;
   }else{
      return 0;
   }
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
my ($file,$line)=@_;
#print "$file\n";
open OUT, ">$file" or die "$!";
     print OUT "$line";
close OUT;
}

sub parsealign
{
my ($file)=@_;
#### Identity:       7/9 (77.8%)
open FD, "$file" or die "$!";
while(<FD>){
   chomp $_;
   my $line=$_;
   if ($line=~/# Identity:\s+(\d+)\/(\d+)\s*/){
      return ($1,$2);
   }
}
close FD;

}

sub align1
{
my ($seq1,$seq2)=@_;
my $test1=">test1\n$seq1\n";
my $test2=">test2\n$seq2\n";
my $file1="test1.txt";
my $file2="test2.txt";
writefile($file1,$test1);
writefile($file2,$test2);
`/home/biosoftware/EMBOSS-5.0.0/emboss/water test1.txt test2.txt -gapopen 10 -gapextend 0.5 -brief Y -outfile test.align`;
my ($identity,$length)=parsealign("test.align");
return ($identity,$length);
}


