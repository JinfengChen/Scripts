#!/usr/bin/perl
use Getopt::Long;
use Statistics::R;

GetOptions (\%opt,"table:s","number:s","help");


my $help=<<USAGE;
Do Chi sequare Test for a file.

Example: 
name	#Test in FF	#Ref in FF	#Test in Rice	#Ref in Rice
apoplast	30	20	39	30
carbohydrate binding	60	68	143	138


Output:
name    #Test in FF     #Ref in FF      #Test in Rice   #Ref in Rice    P-value FDR     
apoplast        30      20      39      30      0.7044  1.0000000
carbohydrate binding    60      68      143     138     0.4515  1.0000000

Run: perl ChisqTestR.pl -table fortest.txt -number 2 > pfam.test
-table: table for test
-number: number of object(species) to compare
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 


my $R=Statistics::R->new();
my @obs;
my @p;
my @fdr;
my $head;
open IN, "$opt{table}" or die "$!";
$R->start_sharedR;
while(<IN>){
    chomp $_;
    my $line=$_;
    #print "$_\n";
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[1]=~/\D+/){
        #print "$line\tP-value\tFDR\tTest\n";
        $head="$line\tP-value\tFDR\t";
        next;
    }
    my $NA=0;
    my @data;
    for(my $i=1;$i<=$opt{number}*2;$i++){
       if ($unit[$i] eq ""){
          $NA=1;
       }
       push (@data,$unit[$i]);
    }
    if ($NA==1){
        push (@obs,$line);
        push (@p,"NA");
        #print "$line\tNA\tNA\tNA\n";
        next;
    }
    my $test=join(",",@data);
    $R->send("x <- matrix(c($test),nc=$opt{number})");  
    $R->send("print (x)");
    my $m=$R->read; 
    #print "$m\n";
    $R->send("y <- chisq.test(x,correct=FALSE)");
    $R->send("print (y)");
    my $y=$R->read;
    #print "$y\n";
    my $pvalue;
    if ($y=~/p-value = (.*)$/ or $y=~/p-value < (.*)$/ or $y=~/p-value > (.*)$/){
       $pvalue=$1;
    }    
    #print "$pvalue\n";
    push (@obs,$line);
    push (@p,$pvalue);
}
close IN;
#################Adjust P-values for Multiple Comparisons############
my $temp=join(",",@p);
#print "$temp\n";
$R->send("parray <- c($temp)");
#$R->send("print (parray)");
#my $a=$R->read;
#print "P:$a\n";
$R->send("q <- p.adjust(parray, method='bonferroni')");
$R->send("print (q)");
my $q=$R->read;
#print "Q:$q\n";
#################delete some [1] element in array###################
@fdr=split(" ",$q);
for(my $j=0;$j<@fdr;$j++){
   if ($fdr[$j]=~/\[|\]/){
     splice (@fdr,$j,1);
   }
}
#################output ############################################

print "$head\n";
for(my $i=0;$i<@obs;$i++){
    print "$obs[$i]\t$p[$i]\t$fdr[$i]\n";
}
$R->stopR();

