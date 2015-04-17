#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"project1:s","project2:s","help");


my $help=<<USAGE;
perl $0 --project1 --project2
--project1,2: file prefix for gene statistic, OBa for OBa.share.cds_length.txt.

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

compareStructureDrawLine($opt{project1},$opt{project2});

sub compareStructureDrawLine
{
my ($title1,$title2)=@_;
open OUT, ">$title1.VS.$title2.r" or die "$!";
print OUT <<"END.";
read.table("$title1.share.cds_length.txt") -> x1
read.table("$title1.share.exon_length.txt") -> x2
read.table("$title1.share.exon_number.txt") -> x3
read.table("$title1.share.intron_length.txt") -> x4
read.table("$title1.share.mRNA_length.txt") -> x5
read.table("$title1.share.intergenic_length.txt") -> x6

read.table("$title2.share.cds_length.txt") -> y1
read.table("$title2.share.exon_length.txt") -> y2
read.table("$title2.share.exon_number.txt") -> y3
read.table("$title2.share.intron_length.txt") -> y4
read.table("$title2.share.mRNA_length.txt") -> y5
read.table("$title2.share.intergenic_length.txt") -> y6


plotnumberline <- function(step,x1,y1,xlim,title,xlab,ylab){
max <- max(x1,y1)+step
br <- seq(0,max,by=step)
tx <-hist(x1,breaks=br,plot=F)
ty <-hist(y1,breaks=br,plot=F)
ylim <- max(tx\$count,ty\$count)+500
xvalue=tx\$mid+step/2
plot(xvalue,tx\$count,ylim=c(0,ylim),xlim=c(0,xlim),main=title,xlab=xlab,ylab=ylab,type="b",pch=1,col=2)
lines(xvalue,ty\$count,ylim=c(0,ylim),xlim=c(0,xlim),type="b",pch=2,col=3)
legend("topright",c("$title1","$title2"),lty=1,pch=1:2,lwd=1,col=c(2,3))
}


plotlengthline <- function(step,x1,y1,xlim,title,xlab,ylab){
max <- max(x1,y1)+step
br <- seq(0,max,by=step)
tx <-hist(x1,breaks=br,plot=F)
ty <-hist(y1,breaks=br,plot=F)
ylim <- max(tx\$count,ty\$count)+500
xvalue=tx\$mid
plot(xvalue,tx\$count,ylim=c(0,ylim),xlim=c(0,xlim),main=title,xlab=xlab,ylab=ylab,type="b",pch=1,col=2)
lines(xvalue,ty\$count,ylim=c(0,ylim),xlim=c(0,xlim),type="b",pch=2,col=3)
legend("topright",c("$title1","$title2"),lty=1,pch=1:2,lwd=1,col=c(2,3))
}

pdf("$title1.VS.$title2.pdf")

plotlengthline(200,x1[,2],y1[,2],5000,"cds_length","Length (bp)","Count")
plotlengthline(200,x2[,2],y2[,2],1000,"exon_length","Length (bp)","Count")
plotnumberline(1,x3[,2],y3[,2],10,"exon_number","Exon Number","Count")
plotlengthline(200,x4[,2],y4[,2],1000,"intron_length","Length (bp)","Count")
plotlengthline(200,x5[,2],y5[,2],5000,"mRNA_length","Length (bp)","Count")
plotlengthline(5000,x6[,2],y6[,2],100000,"intergenic_length","Length (bp)","Count")
dev.off()
END.
close OUT;
system ("cat $title1.VS.$title2.r | R --vanilla --slave");
}
 
