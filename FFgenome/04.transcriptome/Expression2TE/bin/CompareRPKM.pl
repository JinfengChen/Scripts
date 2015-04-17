#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","rpkm1:s","rpkm2:s","project:s","help");


my $help=<<USAGE;
Compare RPKM from different experiments or tissue for one species.
Drow dot plot of RPKM of different tissue for one species.
-bed: BED format gene annotation file, generate by GFF2BED.pl from GFF3 format gene annotation.
-rpkm1: gene expression data 1.
-rpkm2: gene expression data 2.
Gene	RPKM
OBR_GLEAN_10017382 64.3014830721204
OBR_GLEAN_10011719 78.1501432555481
OBR_GLEAN_10007059 24.7132214455825
-project: project name.

Run: perl $0 -bed FF.mRNA.bed -rpkm1 ../input/FF.shoot.rpkm.bgi -rpkm2 ../input/FF.shoot.rpkm.tophat -project BGI_tophat

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $BEDtools="/home/jfchen/software/tools/BEDTools/bin";
my $refrpkm1=rpkmexpr($opt{rpkm1});
my $refrpkm2=rpkmexpr($opt{rpkm2});

open OUT, ">$opt{project}.4r" or die "$!";
open IN, "$opt{bed}" or die "$!"; 
    while(<IN>){
      chomp $_;
      next if ($_ eq "");
      my @unit=split("\t",$_); 
      if (exists $refrpkm1->{$unit[3]} and exists $refrpkm2->{$unit[3]}){
         print OUT "$unit[3]\t$refrpkm1->{$unit[3]}\t$refrpkm2->{$unit[3]}\n";
      }elsif(exists $refrpkm1->{$unit[3]}){
         print OUT "$unit[3]\t$refrpkm1->{$unit[3]}\tNA\n";
      }elsif(exists $refrpkm2->{$unit[3]}){
         print OUT "$unit[3]\tNA\t$refrpkm2->{$unit[3]}\n";
      }else{
         print OUT "$unit[3]\tNA\tNA\n";
      }
    }

close IN;
close OUT;

my @axes=split("\_",$opt{project});

open OUT, ">$opt{project}.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.4r") -> x
pdf("$opt{project}.pdf")
plot(x[,2],x[,3],xlab="$axes[0] RPKM",ylab="$axes[1] RPKM",xlim=c(0,2000),ylim=c(0,2000))
plot(log2(x[,2]),log2(x[,3]),xlab="$axes[0] log2(RPKM)",ylab="$axes[1] log2(RPKM)",xlim=c(0,20),ylim=c(0,20))
hist(log(x[,2]),xlim=c(-10,20),breaks=100,xlab="$axes[0] log(RPKM)",col=2)
hist(log2(x[,2]),xlim=c(-10,20),breaks=100,xlab="$axes[0] log2(RPKM)",col=2)
hist(log10(x[,2]),xlim=c(-10,20),breaks=100,xlab="$axes[0] log10(RPKM)",col=2)
dev.off()
END.
close OUT;

system ("cat $opt{project}.r | R --vanilla --slave");
################################################


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
       #$hash{$unit[0]}=0;
       next;
    }else{
       $hash{$unit[0]}=$unit[1];
    }
}
close IN;
return \%hash;
}




 
