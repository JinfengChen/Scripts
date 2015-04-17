#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","gff:s","rpkm:s","windows:s","project:s","help");


my $help=<<USAGE;
Relation of TE density and Expression levels for gene in one species.
-bed: BED format gene annotation file, generate by GFF2BED.pl from GFF3 format gene annotation.
-gff: BED format TE annotation file, should be contain only TE class that interested to calculate correlation.
-rpkm: gene expression data.
Gene	RPKM
OBR_GLEAN_10017382	64.3014830721204
OBR_GLEAN_10011719	78.1501432555481
OBR_GLEAN_10007059	24.7132214455825
-windows: windows size upstream/downstream of gene, in which we calculate the fraction of TE.
-project: project name.

Run: perl Expression2TE.pl -bed FF.mRNA.bed -gff ../input/chr01/FF.COPIA.bed -rpkm ../input/FF.shoot.rpkm -windows 50000 -project A01_Gene_COPIA
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
my $refrpkm=rpkmexpr($opt{rpkm});
system("$BEDtools/windowBed -a $opt{bed} -b $opt{gff} -w $opt{windows} > $opt{project}.windowsBED");
my $refBED=windowsBED("$opt{project}.windowsBED");
`rm *.windowsBED`;

#& dotplot($refBED,$refrpkm);
& barplot($refBED,$refrpkm);

sub barplot
{
my ($refBED,$refrpkm)=@_;
my %hash;
open OUT, ">$opt{project}.gene_rpkm_TE.density" or die "$!";
foreach(keys %$refBED){
    my $freq=$refBED->{$_}/(2*$opt{windows});
    $freq = $freq > 1 ? 1 : $freq;
    my $rpkm= $refrpkm->{$_} ? $refrpkm->{$_} : "NA";
    print OUT "$_\t$rpkm\t$freq\n";
    my $index=int ($freq*10);
    next if $index > 8;
    #print "Index\t$index\tRPKM\t$rpkm\n";
    push (@{$hash{$index}},$rpkm);
}
close OUT;

my $kb=2*$opt{windows}/1000;

open OUT, ">$opt{project}.bar.4r" or die "$!";
foreach(sort {$a <=> $b} keys %hash){
    my $percent=$_*10;
    my ($mean,$se,$number)=ci($hash{$_});
    my $cilow=$mean-$se;
    my $cihigh=$mean+$se;
    print OUT "$mean\t$cilow\t$cihigh\t$percent\t$se\n";
}
close OUT;

open OUT, ">$opt{project}.bar.r" or die "$!";
print OUT <<"END.";
read.table("$opt{project}.bar.4r") -> x
pdf("$opt{project}.bar.pdf")

#library("gplots")
#barplot2(x[,1],ylim=c(0,60),plot.ci=TRUE,ci.u=x[,3],ci.l = x[,2],names.arg=x[,4],ylab="Expression level (RPKM)",xlab="TE/$kb kb (%)")
#abline(h=0)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-upper, angle=90, code=3, length=length, ...)
}

barx <- barplot(x[,1],ylim=c(-0.5,0.5),names.arg=x[,4],ylab="Normalized RPKM",xlab="TE/$kb kb (%)")
error.bar(barx,x[,1],x[,5]) 
abline(h=0)
dev.off()

END.
close OUT;

system ("cat $opt{project}.bar.r | R --vanilla --slave");
}



sub dotplot
{
my ($refBED,$refrpkm)=@_;
open OUT, ">$opt{project}.dot.4r" or die "$!";
foreach(keys %$refBED){
    my $freq=$refBED->{$_}/(2*$opt{windows});
    $freq = $freq > 1 ? 1 : $freq;
    #print "$_\t$refrpkm->{$_}\n";
    if (exists $refrpkm->{$_} and $refrpkm->{$_} > 3){
       print OUT "$_\t$refrpkm->{$_}\t$freq\n";
    }else{
       print OUT "$_\tNA\t$freq\n";
    }

}
close OUT;

my $kb=2*$opt{windows}/1000;

open OUT, ">$opt{project}.dot.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.dot.4r") -> x
z <- data.frame(x[,2],x[,3])
pdf("$opt{project}.dot.pdf")
plot(x[,2]~x[,3],ylim=c(0,400),ylab="Expression (RPKM)",xlab="TE/$kb kb",data=z)
abline(lm(x[,2]~x[,3],data=z))
dev.off()
END.
close OUT;

system ("cat $opt{project}.dot.r | R --vanilla --slave");
}

################################################
sub windowsBED
{
#### sum len of feature that fill in windows of gene
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $len=$unit[10]-$unit[9]+1;
    if (exists $hash{$unit[3]}){
        $hash{$unit[3]}+=$len;
    }else{
        $hash{$unit[3]} =$len;
    }
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
    $hash{$unit[0]}=$unit[2];
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
    if ($unit[1]=~/\-/){
       $hash{$unit[0]}="NA";
    }else{
       $hash{$unit[0]}=$unit[1];
    }
}
close IN;
return \%hash;
}
=cut

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
return ($mean,$se,$number);
}


 
