#!/usr/bin/perl

## perl simulateSeq.pl filename insertsize coverage
## perl simulateSeq.pl chr04.txt 200 40
if (@ARGV < 3){
  print "Get size+_100 paired sequences\n";
  print "Usage: perl simulateSeq.pl filename insertsize coverag"; 
  print "perl simulateSeq.pl chr04.txt 200 40";
  exit;
}


my $file=$ARGV[0];
my $insertsize=$ARGV[1];
my $coverage=$ARGV[2];
$/=">";
print "$file\n";
open FILE, "$file" or die "can not open my input file";
while (<FILE>) {
	 my @unit=split("\n",$_);
	 my $headline=shift @unit;
	 my $clone;
	 if ($headline=~/(\w+)/) {
            $clone=$1;
	 }
	 my $seq=join("",@unit);
         if ($clone=~/(\w+)/) {
	    randomPEreads ($seq,$clone,$insertsize,$coverage);
	 }
	 my $clonelength=length $seq;
     #print "$headline\n$clonelength\n";
}
close FILE;

sub randomPEreads{

my $mer=75;
my ($rawseq,$clone,$insize,$covarage)=@_;
#print "$clone\n$rawseq\n";
my $revseq=reverse $rawseq;
$revseq=~tr/atgcATGC/tacgTACG/;
my $clonelength=length $rawseq;
print "$clonelength\n";
my $varation=0.2*$insize; ## total varation of insert size, including shorter and larger
my $numberofread=$clonelength*$coverage/(2*$mer);
my $usefullseq=$clonelength-$insize-50-2*$mer+1; ## used to generate position for first read
my $counter;                                     ## need to ensure the second read is in the sequence
my $outfile=$clone."_".$insize."_".$coverage."_".$mer."_".reads1.".".fa;
my $outfile2=$clone."_".$insize."_".$coverage."_".$mer."_".reads2.".".fa;
print "$outfile\t$outfile2\n";
open READ, ">$outfile" or die "can not open my outfile read";
open READ2, ">$outfile2" or die "can not open my outfile read";
for (my $i=0;$i<$numberofread;$i++) {
    my $strand=rand 1;
	if ($strand <= 0.5) {
	   $counter++;
	   my $index=int (rand $usefullseq);
	   #print "$index\n";
	   my $insertion=int (rand $varation)+0.9*$insize-75;
	   my $index2=$index+$insertion;
	   my $read=substr($rawseq,$index,$mer);
	   my $read2=substr($rawseq,$index2,$mer);
           if ($insize < 2000){
              $read2=reverse $read2;
              $read2=~tr/atgcATGC/tacgTACG/;
           }else{
              $read=reverse $read;
              $read=~tr/atgcATGC/tacgTACG/;
           }
	   print READ ">READ\_$counter\_f\n$read\n";
           print READ2 ">READ\_$counter\_r\n$read2\n";

    }else{
	   $counter++;
	   my $index=int (rand $usefullseq);
	   #print "$index\n";
           #my $insertion=int (rand 200)+$insize-100-75;
           my $insertion=int (rand $varation)+0.9*$insize-75;
           my $index2=$index+$insertion;
	   my $read=substr($revseq,$index,$mer);
           my $read2=substr($revseq,$index2,$mer);
           if ($insize < 2000){
               $read2=reverse $read2;
               $read2=~tr/atgcATGC/tacgTACG/;
           }else{
              $read=reverse $read;
              $read=~tr/atgcATGC/tacgTACG/;
           }
	   print READ2 ">READ\_$counter\_r\n$read\n";
           print READ ">READ\_$counter\_f\n$read2\n";
	
	}
}
close READ;
close READ2;
} 

sub randomreads{

my $coverage=10;
my $mer=75;
my ($rawseq,$clone)=@_;
#print "$clone\n$rawseq\n";
my $revseq=reverse $rawseq;
$revseq=~tr/atgcATGC/tacgTACG/;
my $clonelength=length $rawseq;
print "$clonelength\n";
my $numberofread=$clonelength*$coverage/$mer;
my $usefullseq=$clonelength-$mer+1;
my $counter;
my $outfile=$clone."_".$coverage."_".$mer."_".reads.".".fa;
print "$outfile\n";
open READ, ">$outfile" or die "can not open my outfile read";
for (my $i=0;$i<$numberofread;$i++) {
    my $strand=rand 1;
	if ($strand <= 0.5) {
	   $counter++;
	   my $index=int (rand $usefullseq);
	   #print "$index\n";
	   my $read=substr($rawseq,$index,$mer);
	   print READ ">READ\_$counter\_f\n$read\n";
    }else{
	   $counter++;
	   my $index=int (rand $usefullseq);
	   #print "$index\n";
	   my $read=substr($revseq,$index,$mer);
	   print READ ">READ\_$counter\_r\n$read\n";
	
	}
}
close READ;
} 
