#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"ortholog:s","rpkm1:s","rpkm2:s","project:s","help");


my $help=<<USAGE;
Compare RPKM of ortholog genes bewteen two species.
Draw dot plot for expression levels of ortholog genes between two species FF and rice.
-ortholog: ortholog table of two species
-rpkm1: gene expression data 1.
-rpkm2: gene expression data 2.
Gene	RPKM
OBR_GLEAN_10017382 64.3014830721204
OBR_GLEAN_10011719 78.1501432555481
OBR_GLEAN_10007059 24.7132214455825
-project: project name.

Run: perl $0 -ortholog OB-OS.orth -rpkm1 ../input/FF.shoot.rpkm.tophat -rpkm2 ../input/RAP3.shoot.rpkm.tophat -project FF_Rice 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
my $refrpkm1=rpkmexpr($opt{rpkm1});
my $refrpkm2=rpkmexpr($opt{rpkm2});

open OUT, ">$opt{project}.4r" or die "$!";
open IN, "$opt{ortholog}" or die "$!"; 
    while(<IN>){
      chomp $_;
      next if ($_ eq "");
      my @unit=split("\t",$_); 
      if (exists $refrpkm1->{$unit[0]} and exists $refrpkm2->{$unit[1]}){
         print OUT "$unit[0]\t$refrpkm1->{$unit[0]}\t$refrpkm2->{$unit[1]}\n";
      }elsif(exists $refrpkm1->{$unit[0]}){
         print OUT "$unit[0]\t$refrpkm1->{$unit[0]}\tNA\n";
      }elsif(exists $refrpkm2->{$unit[3]}){
         print OUT "$unit[0]\tNA\t$refrpkm2->{$unit[1]}\n";
      }else{
         print OUT "$unit[0]\tNA\tNA\n";
      }
    }

close IN;
close OUT;

my @axes=split("\_",$opt{project});

open OUT, ">$opt{project}.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.4r") -> x
pdf("$opt{project}.pdf")
plot(x[,2],x[,3],xlab="$axes[0] RPKM",ylab="$axes[1] RPKM")
hist(x[,2],xlim=c(-5,5),breaks=100,xlab="$axes[0] RPKM",col=2)
hist(x[,3],xlim=c(-5,5),breaks=100,xlab="$axes[1] RPKM",col=3)
hist(x[,2],xlim=c(-5,5),breaks=100,xlab="$axes[0] RPKM",col=2)
par(new=TRUE)
hist(x[,3],xlim=c(-5,5),breaks=100,xlab="",ylab="",main="",axes=FALSE,new=T,col=3)
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
    $hash{$unit[0]}=$unit[2];
}
close IN;
return \%hash;
}




 
