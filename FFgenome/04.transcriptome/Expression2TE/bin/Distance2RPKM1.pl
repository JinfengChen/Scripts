#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","gff:s","rpkm:s","project:s","help");


my $help=<<USAGE;
Relation of distance from nearest TE and Expression level for gene in one species.
-bed: BED format gene annotation file, generate by GFF2BED.pl from GFF3 format gene annotation.
-gff: BED format TE annotation file, should be contain only TE class that interested to calculate correlation.
-rpkm: gene expression data.
Gene	RPKM
OBR_GLEAN_10017382	64.3014830721204
OBR_GLEAN_10011719	78.1501432555481
OBR_GLEAN_10007059	24.7132214455825
-project: project name.

Run: perl $0 -bed FF.mRNA.bed -gff ../input/FF.repeat.gff -rpkm ../input/FF.shoot.rpkm -project FF
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash; ## change this hash to meth and unmeth to store distance and expression
my @control;
my $control;
my $bin=200;
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
my $refrpkm=rpkmexpr($opt{rpkm});
system("$BEDtools/closestBed -a $opt{bed} -b $opt{gff} > $opt{project}.closestBED");
my $refBED=closestBED("$opt{project}.closestBED");
`rm *.closestBED`;
open OUT, ">$opt{project}.4r" or die "$!"; ## open two file one for meth and the other for unmeth
foreach(keys %$refBED){
    $rpkm= $refrpkm->{$_} ? $refrpkm->{$_} : "NA";
    next if ($rpkm eq "NA");
    my $index;
    if ($refBED->{$_} > 0){
       $index=int ($refBED->{$_}/$bin) + 1;
    }else{
       $index=0;
    }
    my $dist=$index*$bin;
    next if ($dist > 2500);
    push (@{$hash{$index}},$rpkm); ## test if methylation or not and push to meth or unmeth
    if ($index > 1){
       push (@control,$rpkm);
    }
}
$control=mean(\@control);
foreach(sort {$a <=> $b} keys %hash){
    my $array=$hash{$_};
    my $dist =$_*$bin;
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;


open OUT, ">$opt{project}.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.4r") -> x
pdf("$opt{project}.pdf")
library("gplots")
plotmeans(x[,2]~x[,1],ylim=c(0.5,1.5),connect=FALSE,ylab="Expression (RPKM)",xlab="Distance to gene (bp)")
abline(h=$control)
dev.off()
END.
close OUT;

system ("cat $opt{project}.r | R --vanilla --slave");
################################################
sub closestBED
{
#### get the distance from gene to nearest TE 
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $distance;
    if ($unit[9] > $unit[2]){
       $distance=$unit[9]-$unit[2]+1;
    }elsif($unit[9] < $unit[1]){
       $distance=$unit[1]-$unit[9]+1;
    }else{
       $distance=0;
    }
    $hash{$unit[3]}=$distance;
}
close IN;
return \%hash;
}



sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split(" ",$_);
    if ($unit[1]=~/\-/ or $unit[1] == 0){
       $hash{$unit[0]}="NA";
    }else{
       $hash{$unit[0]}=log10($unit[1]);
    }
}
close IN;
return \%hash;
}

sub mean
{
my ($num)=@_;
my $loop=0;
my $total;
foreach  (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
$number=@$num;
$mean=$total/$loop;
$median=$num->[int $number/2];
return $mean;
}

sub log10 {
    my ($n) = shift;
    return log($n)/log(10);
}
 
