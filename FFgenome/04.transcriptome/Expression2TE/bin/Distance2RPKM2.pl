#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","gff:s","rpkm:s","me:s","project:s","help");


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
-me: methylation status(id,status)
Run: perl $0 -bed FF.mRNA.bed -gff ../input/FF.repeat.gff -rpkm ../input/FF.shoot.rpkm -me ../input/FF.repeat.gff.Me.status -project FF
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash; ## change this hash to meth and unmeth to store distance and expression
my %me;
my %unme;
my @control;
my $control;
my $max=1000;
my $bin=$max/10;
my $controlindex=500/$bin;
my $cutoff=0.1;
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
my $refrpkm=rpkmexpr($opt{rpkm});
my $refme  =mestatus($opt{me});
system("$BEDtools/closestBed -a $opt{bed} -b $opt{gff} > $opt{project}.closestBED");
my $refBED=closestBED("$opt{project}.closestBED");
`rm *.closestBED`;
open OUT, ">$opt{project}.gene.closestTE.inf" or die "$!";
print OUT "Gene\tRPKM\tRepeat\tDistance\tMethylation\n";
foreach(keys %$refBED){
    my $rpkm0= $refrpkm->{$_}->[0] ? $refrpkm->{$_}->[0] : "NA";
    my $rpkm= $refrpkm->{$_}->[1] ? $refrpkm->{$_}->[1] : "NA";
    print OUT "$_\t$rpkm0\t$refBED->{$_}->[0]\t$refBED->{$_}->[1]\t$refme->{$refBED->{$_}->[0]}\n";
    next if ($rpkm eq "NA");
    my $index;
    if ($refBED->{$_}->[1] > 0){
       $index=int ($refBED->{$_}->[1]/$bin) + 1;
    }else{
       $index=0;
    }
    my $dist=$index*$bin;
    next if ($dist > $max);
    push (@{$hash{$index}},$rpkm); ## test if methylation or not and push to meth or unmeth
    if ($refme->{$refBED->{$_}->[0]} >= $cutoff){
         push (@{$me{$index}},$rpkm);
    }elsif($refme->{$refBED->{$_}->[0]} < $cutoff and $refme->{$refBED->{$_}->[0]} ne "NA"){
         push (@{$unme{$index}},$rpkm);
    }
    if ($index > $controlindex){
       push (@control,$rpkm);
    }
}
close OUT;
$control=mean(\@control);

open OUT, ">$opt{project}.all.4r" or die "$!";
my @pos;
my @all;
foreach(sort {$a <=> $b} keys %hash){
    my $array=$hash{$_};
    my $temp=@$array;
    my $dis =$_*$bin;
    push (@pos,$dis);
    push (@all,$temp);
    my $dist =$_*$bin;
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;

my @me;
open OUT, ">$opt{project}.me.4r" or die "$!";
foreach(sort {$a <=> $b} keys %me){
    my $array=$me{$_};
    my $dist =$_*$bin;
    my $temp=@$array;
    push (@me,$temp);
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;

my @unme;
open OUT, ">$opt{project}.unme.4r" or die "$!";
foreach(sort {$a <=> $b} keys %unme){
    my $array=$unme{$_};
    my $temp=@$array;
    push (@unme,$temp);
    my $dist =$_*$bin;
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;


open SUM, ">$opt{project}.summary" or die "$!";
print SUM "Distance\tAll\tMe\tUnMe\n";
for(my $i=0;$i<@all;$i++){
   print SUM "$pos[$i]\t$all[$i]\t$me[$i]\t$unme[$i]\n";
}
close SUM;


open OUT, ">$opt{project}.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.all.4r") -> x
read.table("$opt{project}.me.4r") -> y
read.table("$opt{project}.unme.4r") -> z
pdf("$opt{project}.pdf")
library("gplots")
plotmeans(y[,2]~y[,1],ylim=c(-1,1),n.label=FALSE,barcol="black",ylab="Expression (RPKM)",xlab="Distance to gene (bp)")
par(new=T)
plotmeans(z[,2]~z[,1],ylim=c(-1,1),n.label=FALSE,col="red",barcol="red",xlab="",ylab="",new=TRUE);
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
    my $id;
    if ($unit[9] > $unit[2]){
       $distance=$unit[9]-$unit[2]+1;
    }elsif($unit[10] < $unit[1]){
       $distance=$unit[1]-$unit[10]+1;
    }else{
       $distance=0;
    }
    if ($unit[14]=~/ID=(.*?);/){
       $id=$1;
    }
    $hash{$unit[3]}=[$id,$distance];
}
close IN;
return \%hash;
}

sub mestatus
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

sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    $hash{$unit[0]}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
}

=pod
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
       $hash{$unit[0]}=[-2,0.01];
       #$hash{$unit[0]}=["NA","NA"];
    }else{
       my $temp=log10($unit[1]);
       $hash{$unit[0]}=[$temp,$unit[1]];
    }
}
close IN;
return \%hash;
}
=cut

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
 
