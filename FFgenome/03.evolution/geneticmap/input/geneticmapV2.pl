#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"map:s","help");


my $help=<<USAGE;
perl $0 --map ../input/MSU.genetic.map
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $win=1000000;
my $step=500000;

###read map and store in %map;
my %map;
my %hash;
open IN, "$opt{map}" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   next if ($unit[7] eq "NS");
   my $chr="Chr$unit[3]";
   my $cm;
   if ($unit[6]=~/$chr\;(.*)cM/){
      $cm=$1;
   }else{
      next;
   }
   push (@{$map{$unit[3]}},[$unit[4],$cm]);
}
close IN;

###filter out duplicated and wrong markers
foreach (sort {$a <=> $b} keys %map){
    my $dist=$map{$_};
    print "Chr$_\n";
    my $chr="Chr$_";
    my %temp;
    foreach(@$dist){
       $temp{$_->[0]}=$_->[1];
    }
    my $last=0;
    my %final;
    foreach(sort {$a <=> $b} keys %temp){
        next if ($temp{$_} < $last);
        $last=$temp{$_};
        #print "$_\t$temp{$_}\n";
        $final{$_}=$temp{$_};
    }
    #####
    my $title="$chr RecombinationRate";
    my $Rfile="$chr\.RecombinationRate.4r";
    my $Rdraw="$chr\.RecombinationRate.r";
    my $Rpdf="$chr\.RecombinationRate.pdf";
    open OUT, ">$Rfile" or die "$!";
    my %wininf;
    foreach(sort {$a <=> $b} keys %final){
       my $line="$_\t$final{$_}";
       #print "$line\n";
       sliding($line,$win,$step,\%wininf);
    }
    
    foreach(sort {$a <=> $b} keys %wininf){
        print "$_\n";
        my @sortinf=sort {$a->[0] <=> $b->[0]} @{$wininf{$_}};
        foreach(@sortinf){
           print "$_->[0]\t$_->[1]\n";
        }
        if (@sortinf > 1 ){
           my $rr=1000000*($sortinf[$#sortinf][1]-$sortinf[0][1])/($sortinf[$#sortinf][0]-$sortinf[0][0]);
           print OUT "$_\t$rr\n";
        }else{
           print OUT "$_\t0\n";
        }
    }
    close OUT; 
    & drawrr($title,$Rfile,$Rdraw,$Rpdf);
}

####

sub sliding
{
my ($line,$win,$step,$hash)=@_;
my @unit=split("\t",$line);
my $bin=int ($unit[0]/$step) + 1;
for(my $i=$bin;$i>0;$i--){
   my $start=($i-1)*$step;
   my $end=$start+$win;
   #print "$i\t$unit[0]\t$start\t$end\n";
   if ($unit[0] >= $start and $unit[0] < $end) {
            push (@{$hash->{$start}},[$unit[0],$unit[1]]);
            #print "$i\t$unit[0]\t$start\t$end\n";
   }elsif($unit[0] >= $end or $unit[0] < $start){
            return;
   }
}
}


sub drawrr
{
my ($title,$file,$draw,$pdf)=@_;
open OUT, ">$draw" or die "$!";
print OUT <<"END.";
read.table("$file") -> x
pdf("$pdf")
barplot(x[,2],xlab="Position On Chromosome",ylab="Recombination Rate (cM/Mb)",main="$title");
plot(x[,1],x[,2],type="l",xlab="Position On Chromosome",ylab="Recombination Rate (cM/Mb)",main="$title");
dev.off()
END.
close OUT;
system ("cat $draw | R --vanilla --slave");
}



