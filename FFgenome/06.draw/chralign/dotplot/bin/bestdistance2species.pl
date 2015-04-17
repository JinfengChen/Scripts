#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"head:s","help");


my $help=<<USAGE;
summary best distance of non-syntenic duplicated gene.
perl $0 --head
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable("$opt{head}.os.bestdistance");
readtable("$opt{head}.ob.bestdistance");
draw2compare("$opt{head}");

sub readtable
{
my ($file)=@_;
my $head=$1 if ($file=~/(.*)\.bestdistance/);
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if ($unit[2] <= 10){
       $hash{"10"}++;
    }elsif($unit[2] <= 100){
       $hash{"100"}++;
    }elsif($unit[2] <= 1000){
       $hash{"1000"}++;
    #}elsif($unit[2] == 10000){
    #   $hash{"10000"}++;
    }else{
       $hash{"5000"}++;
    }
}
close IN;
my %note=(
   "10" => "<10",
   "100"=> "10-100",
   "1000"=>"100-1000",
   #"5000"=>"1000-5000",
   "5000"=>">1000"
);
open OUT, ">$head.bestdistance.4r" or die "$!";
foreach (sort {$a <=> $b} keys %hash){
    print OUT "$_\t$note{$_}\t$hash{$_}\n";
}
close OUT;
}

sub draw2compare
{
my ($head)=@_;
open OUT, ">$head.bestdistance.r" or die "$!";
print OUT <<"END.";
read.table("$head.os.bestdistance.4r") -> x
read.table("$head.ob.bestdistance.4r") -> y
data=rbind(x[,3],y[,3])
pdf("$head.bestdistance.pdf")
barplot(data,beside=TRUE,border=NA,names.arg=x[,2],xlab="Distance to best hit",ylab="Gene Number",col=3:4);
par(font=3)
legend("topleft",legend=c("Oryza sativa","Oryza brachyantha"),fill=c(3,4),bty="n")
dev.off()
END.
close OUT;
system ("cat $head.bestdistance.r | R --vanilla --slave");

}
 
