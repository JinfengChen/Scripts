#!/usr/bin/perl
use Getopt::Long;
use Statistics::R;

GetOptions (\%opt,"table:s","help");


my $help=<<USAGE;
Do Chi sequare Test for a file.
Example: 
name	#Test in FF	#Ref in FF	#Test in Rice	#Ref in Rice
apoplast	30	20	39	30
carbohydrate binding	60	68	143	138

T is the expected number for each cell
n is the total number of 4 cell
If have a T < 5 or n < 40, we use fisher exactly test. else we use Chi Square Test.

Output:
name    #Test in FF     #Ref in FF      #Test in Rice   #Ref in Rice    P-value FDR     
apoplast        30      20      39      30      0.7044  1.0000000
carbohydrate binding    60      68      143     138     0.4515  1.0000000

Run: perl ChisqTestR.pl -table ../input/fortest.txt > fortest.log
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
open IN, "$opt{table}" or die "$!";;
$R->startR;
while(<IN>){
    chomp $_;
    my $line=$_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[1]=~/\D+/){
        #print "$line\tP-value\tFDR\tTest\n";
        $head="$line\tP-value\tFDR";
        next;
    }
    if ($unit[1] eq "" or $unit[2] eq "" or $unit[3] eq "" or $unit[4] eq ""){
        push (@obs,$line);
        push (@p,"NA");
        #print "$line\tNA\tNA\tNA\n";
        next;
    }
    my $n1_ =$unit[1]+$unit[2];
    my $n2_ =$unit[3]+$unit[4];
    my $n_1 =$unit[1]+$unit[3];
    my $n_2 =$unit[2]+$unit[4];
    my $n__ =$n1_ + $n2_;
    $R->send("x <- matrix(c($unit[1],$unit[3],$unit[2],$unit[4]),nc=2)");
    unless (($n1_*$n_1)/$n__ < 5 or ($n1_*$n_2)/$n__ < 5 or ($n2_*$n_1)/$n__ < 5 or ($n2_*$n_2)/$n__ < 5 or $n__ < 40){ 
       #$R->send("x <- matrix(c($unit[1],$unit[3],$unit[2],$unit[4]),nc=2)");   
       $R->send("y <- chisq.test(x,correct=FALSE)");
       $R->send("print (y)");
       my $y=$R->read;
       #print "$y\n";
       my $pvalue;
       if ($y=~/p-value = (.*)$/ or $y=~/p-value > (.*)$/ or $y=~/p-value < (.*)$/){
          $pvalue=$1;
       }    
       #print "$line\t$pvalue\tChisq\n";
       push (@obs,$line);
       push (@p,$pvalue);
    }else{
       #$R->send("x <- matrix(c($unit[1],$unit[3],$unit[2],$unit[4]),nc=2)");
       $R->send("y <- fisher.test(x)");
       $R->send("print (y)");
       my $y=$R->read;
       my $pvalue;
       if ($y=~/p-value = (.*)\n/ or $y=~/p-value < (.*)\n/ or $y=~/p-value > (.*)\n/){
          $pvalue=$1;
       }
       #print "$line\t$pvalue\tFisher\n";
       push (@obs,$line);
       push (@p,$pvalue);
    }
}
close IN;
#################output ############################################
#print "$head\n";
#for(my $i=0;$i<@obs;$i++){
#    print "$obs[$i]\t$p[$i]\n";
#}
#$R->stopR();


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


