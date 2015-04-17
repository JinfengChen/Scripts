#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;
use warnings;
use strict;

our %opt;

GetOptions(\%opt,"window:s","chrlen:s","bed:s","project:s","help");

$opt{window} ||= 50000;


my $help=<<USAGE;
The BED often contains gap that due to unphase data. We follow the HMMSeg's suggestion to fill these gaps using R.
For gaps less than 2000bp, data are linearly interpolated using the two immediately flanking data points. For gaps larger than 2000bp, we use an adaptive loess fitting strategy. In this case, a linear loess fit is computed at the point of interpolation, using all points in a window of width 50 times the gap to be filled for the fit. The R function loess is used for this purpose, using default weights. The value of the loess fit at the point of interpolation is taken to be the interpolated value there. 
Define chromatin state based on CG,CHG,CHH distribution.
perl $0 --window 50000 --chrlen IRGSP.chrlen -bed rice_methylation_wave -project rice
USAGE

if ($opt{help} or keys %opt < 2){
   print "$help\n";
   exit();
}

my $refchr=chrlen($opt{chrlen});

`rm $opt{bed}/*.fill.*`;
foreach my $chr (keys %$refchr){
    my @file=glob("$opt{bed}/$chr.*.wave.bed");
   `rm $opt{bed}/*local*`;
    foreach my $file (@file){
       print "$file\n";
       my $head;
       if ($file=~/(.*)\.wave\.bed/){
          $head=$1;
       }
       fill($file,$head);
   }
   `rm $opt{bed}/*local*`;
   & list($chr);
   & HMMSeg($chr);

} ##foreach chromosome 
`rm *.list`;
& wig2bar();

#############
sub list
{
my ($chr)=@_;
`find $opt{bed} | grep "CG.wave.fill.bed" | grep $chr > $opt{project}.$chr.cg50000.list`;
`find $opt{bed} | grep "CHG.wave.fill.bed" | grep $chr > $opt{project}.$chr.chg50000.list`;
`find $opt{bed} | grep "CHH.wave.fill.bed" | grep $chr > $opt{project}.$chr.chh50000.list`;
}

sub HMMSeg
{
my ($chr)=@_;
#`java -jar HMMSeg.jar --num-states 4 --input-bed --smooth 100000 --num-starts 10 --log hmm.log $opt{project}.$chr.cg50000.list $opt{project}.$chr.chg50000.list $opt{project}.$chr.chh50000.list > hmm.log1 2> hmm.log2`;
`java -jar HMMSeg.jar --num-states 2 --input-bed --smooth 200000 --num-starts 10 --log hmm.log $opt{project}.$chr.cg50000.list $opt{project}.$chr.chg50000.list $opt{project}.$chr.chh50000.list > hmm.log1 2> hmm.log2`;
}

sub wig2bar
{
`perl wig2bar.pl -wig $opt{bed} > wig.log 2> wig.log2`;
}

############
sub fill
{
my ($file,$head)=@_;
my $out="$head.wave.fill.bed";
open OUT, ">$out" or die "$!";
open IN, "$file" or die "$!";
my $first=<IN>;
chomp $first;
my @last=split("\t",$first);
my @backward;
my @forward;
my %linear;
my %local;
my %pos;
my $index=0;
my @data;
push (@data,[@last]);
while(<IN>){
   chomp $_;
   $index++;
   my $line=$_;
   my @unit=split("\t",$_);
   my $test=($unit[1]-$data[-1][1])/$opt{window};
   if ($test == 1){
   }elsif ($test <= 2){  ### use linear model to fill the data
      #print "Linear\t$test\t$unit[1]\t$data[-1][1]\n$line\n";
      $local{$index}=$test-1;
      $pos{$index}=$data[-1][1];
   }elsif ($test > 2){   ### use local polynomial regression fitting to fill the data
      #print "Losse\t$test\t$unit[1]\t$data[-1][1]\n$line\n";
      $local{$index}=$test-1;
      $pos{$index}=$data[-1][1];
   }
   push (@data,[@unit]);
}

my @add=localfill0(\@data,\%local,$head,\%pos);
push (@data,@add);
@data=sort {$a->[1] <=> $b->[1]} @data;
for(my $i=0;$i<@data;$i++){
   print OUT "$data[$i][0]\t$data[$i][1]\t$data[$i][2]\t$data[$i][3]\t$data[$i][4]\n";
}
close IN;
close OUT;
}

######## fill gap pos with 0
sub localfill0
{
my ($data,$hash,$head,$pos)=@_;
my @add;
my $count;
foreach my $index (keys %$hash){
   my $gap=$hash->{$index};
   my @missing;
   for(my $i=1;$i<=$gap;$i++){
         push (@missing,0);
   }
   my $chr;
   my $type;
   if ($head=~/(chr\d+)\.(\w+)/){
      $chr=$1;
      $type=$2;
   }
   for(my $i=0;$i<@missing;$i++){
       $count++;
       my $name=$type."add".$count;
       my $missstart=$pos->{$index}+($i+1)*$opt{window};
       my $missend  =$missstart+$opt{window}-1;
       push (@add,[$chr,$missstart,$missend,$name,$missing[$i]]);
   }
}
return @add;
}

######## fill gap pos with R model predicted data.
sub localfill
{
my ($data,$hash,$head,$pos)=@_;
my @add;
my $count;
foreach my $index (keys %$hash){
   my $gap=$hash->{$index};
   my $size=25*$gap;
   my $back=getback($data,$index,$size);
   my $forward=getforward($data,$index,$size);
   my $local="$head.$index.local";
   my @missing;
   open TEMP, ">$local" or die "$!";
      my $rank;
      for(my $i=0;$i<@$back;$i++){
         $rank++;
         print TEMP "$rank\t$back->[$i]\n";
      }
      for(my $i=1;$i<=$gap;$i++){
         $rank++;
         print TEMP "$rank\tNA\n";
         push (@missing,$rank);
      }
      for(my $i=0;$i<@$forward;$i++){
         $rank++;
         #print STDERR "$local\n$rank\n";
         print TEMP "$rank\t$forward->[$i]\n";
      }
   close TEMP;
   my $missinglist=join(",",@missing);
   my $rcode=<<CODE;
MissingList <- c($missinglist)
data <- read.table("$local")
pdf("$local.pdf")
x <- data[,1]
y <- data[,2]
plot(x,y)
y.loess <- loess(y ~ x, span=0.75, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
lines(x,y.predict,col=2)
y.Missing <- predict(y.loess, data.frame(x=MissingList))
write.table(y.Missing,"$head.$index.local.R.out");
points(MissingList, y.Missing, pch=FILLED.CIRCLE<-19, col=2)
dev.off()
CODE
   open CODE, ">$head.$index.local.R" or die "$!";
       print CODE "$rcode\n";
   close CODE;
   system("R --vanilla --slave < $head.$index.local.R");
   
   my @missingvalue;
   open MISS, "$head.$index.local.R.out" or die "$!";
       <MISS>;
       while (<MISS>){
           chomp $_;
           my @temp1=split(" ",$_);
           push (@missingvalue,$temp1[1]);
       }
   close MISS;
   
   my $chr;
   my $type;
   if ($head=~/(chr\d+)\.(\w+)/){
      $chr=$1;
      $type=$2;
   } 
   for(my $i=0;$i<@missingvalue;$i++){
       $count++;
       my $name=$type."add".$count;
       my $missstart=$pos->{$index}+($i+1)*$opt{window};
       my $missend  =$missstart+$opt{window}-1;
       push (@add,[$chr,$missstart,$missend,$name,$missingvalue[$i]]);
   }
}
return @add;
}


######
sub getback
{
my ($data,$index,$number)=@_;
my $start=$index-$number < 0 ? 0 : $index-$number;
my @back;
for(my $i=$index-1;$i>=$start;$i--){
   unshift (@back,$data->[$i][4]); 
}
return \@back;
}

######
sub getforward
{
my ($data,$index,$number)=@_;
$index=$index-1;
my $max=@$data-1;
my $end=$index+$number > $max ? $max : $index+$number;
my @forward;
for(my $i=$index+1;$i<=$end;$i++){
   push (@forward,$data->[$i][4]);
}
return \@forward;
}

#################
### get chr length hash
sub chrlen
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}


