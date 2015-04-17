#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"rpkm1:s","rpkm2:s","project:s","help");


my $help=<<USAGE;
Normalization RPKM by standard Z score.
-rpkm1: gene expression data 1.
-rpkm2: gene expression data 2.
Gene	RPKM
OBR_GLEAN_10017382 64.3014830721204
OBR_GLEAN_10011719 78.1501432555481
OBR_GLEAN_10007059 24.7132214455825
-project: project name.

Run: perl $0 -rpkm1 ../input/FF.shoot.rpkm.bgi -rpkm2 ../input/RAP3.shoot.rpkm.tophat -project BGI_tophat

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $refrpkm1=rpkmexpr($opt{rpkm1});
my $refrpkm2=rpkmexpr($opt{rpkm2});

my $normalrpkm1=normalization($refrpkm1);
my $normalrpkm2=normalization($refrpkm2);

open OUT1, ">$opt{project}.1.4r" or die "$!";
foreach(keys %$refrpkm1){
   print OUT1 "$_\t$refrpkm1->{$_}\t$normalrpkm1->{$_}\n";
}
close OUT1;

open OUT2, ">$opt{project}.2.4r" or die "$!";
foreach(keys %$refrpkm2){
   print OUT2 "$_\t$refrpkm2->{$_}\t$normalrpkm2->{$_}\n";
}
close OUT2;

open OUT, ">$opt{project}.r" or die "$!";
print OUT <<"END.";
read.table("$opt{project}.1.4r") -> x1
read.table("$opt{project}.2.4r") -> x2
pdf("$opt{project}.pdf")
hist(x1[,2],breaks=100,xlab="RPKM",col=2)
hist(x1[,3],xlim=c(-5,5),breaks=100,xlab="normalized RPKM",col=2)
hist(x2[,2],breaks=100,xlab="RPKM",col=2)
hist(x2[,3],xlim=c(-5,5),breaks=100,xlab="normalized RPKM",col=2)
hist(x1[,2],breaks=100,xlab="RPKM",col=3)
par(new=TRUE)
hist(x2[,2],breaks=100,xlab="",ylab="",main="",axes=FALSE,col=4,new=T)
hist(x1[,3],xlim=c(-5,5),breaks=100,xlab="RPKM",col=3)
par(new=TRUE)
hist(x2[,3],xlim=c(-5,5),breaks=100,xlab="",ylab="",main="",axes=FALSE,col=4,new=T)
dev.off()
END.
close OUT;

system ("cat $opt{project}.r | R --vanilla --slave");
################################################

sub normalization
{
my ($rpkm)=@_;
my @value=values %$rpkm;
my ($mean,$sd)=ci(\@value);
my %hash;
foreach(keys %$rpkm){
   my $id=$_;
   my $expr=$rpkm->{$id};
   my $zscore=($expr-$mean)/$sd;
   $hash{$id}=$zscore;
}
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
    if ($unit[1]=~/\-/ or $unit[1]==0){
       #$hash{$unit[0]}=-1;
    }else{
       #$hash{$unit[0]}=$unit[1];
       $hash{$unit[0]}=log2($unit[1]);
    }
}
close IN;
return \%hash;
}

sub ci
{
my ($num)=@_;
my $loop=0;
my $total=0;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        #my $temp=log($_);
        #print "$_\t$temp\n";
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}
my $number=$loop;
return (0,0,0) if ($number < 2);
my $mean=$total/$number;
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return ($mean,$SD,$number);
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
return 0 if ($loop == 0);
my $mean=$total/$loop;
return $mean;
}

sub log2 {
    my ($n) = shift;
    return log($n)/log(2);
}



 
