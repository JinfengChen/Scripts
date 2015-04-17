#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"project:s","help");


my $help=<<USAGE;
perl $0 --project
--project: file prefix for gene statistic, OBa for OBa.share.cds_length.txt.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

compareStructureDrawLine($opt{project});

sub compareStructureDrawLine
{
my ($title)=@_;
open OUT, ">$title.r" or die "$!";
print OUT <<"END.";
read.table("$title.share.cds_length.txt") -> x1
read.table("$title.share.exon_length.txt") -> x2
read.table("$title.share.exon_number.txt") -> x3
read.table("$title.share.intron_length.txt") -> x4
read.table("$title.share.mRNA_length.txt") -> x5
read.table("$title.share.intergenic_length.txt") -> x6

read.table("$title.unique.cds_length.txt") -> y1
read.table("$title.unique.exon_length.txt") -> y2
read.table("$title.unique.exon_number.txt") -> y3
read.table("$title.unique.intron_length.txt") -> y4
read.table("$title.unique.mRNA_length.txt") -> y5
read.table("$title.unique.intergenic_length.txt") -> y6

read.table("$title.noncluster.cds_length.txt") -> z1
read.table("$title.noncluster.exon_length.txt") -> z2
read.table("$title.noncluster.exon_number.txt") -> z3
read.table("$title.noncluster.intron_length.txt") -> z4
read.table("$title.noncluster.mRNA_length.txt") -> z5
read.table("$title.noncluster.intergenic_length.txt") -> z6

plotnumberline <- function(step,x1,y1,z1,xlim,title,xlab,ylab){
max <- max(x1,y1,z1)+step
br <- seq(0,max,by=step)
tx <-hist(x1,breaks=br,plot=F)
ty <-hist(y1,breaks=br,plot=F)
tz <-hist(z1,breaks=br,plot=F)
ylim <- max(tx\$count,ty\$count,tz\$count)+500
xvalue=tx\$mid+step/2
plot(xvalue,tx\$count,ylim=c(0,ylim),xlim=c(0,xlim),main=title,xlab=xlab,ylab=ylab,type="b",pch=1,col=2)
lines(xvalue,ty\$count,ylim=c(0,ylim),xlim=c(0,xlim),type="b",pch=2,col=3)
lines(xvalue,tz\$count,ylim=c(0,ylim),xlim=c(0,xlim),type="b",pch=3,col=4)
legend("topright",c("share","unique","noncluster"),lty=1,pch=1:3,lwd=1,col=c(2,3,4))
}


plotlengthline <- function(step,x1,y1,z1,xlim,title,xlab,ylab){
max <- max(x1,y1,z1)+step
br <- seq(0,max,by=step)
tx <-hist(x1,breaks=br,plot=F)
ty <-hist(y1,breaks=br,plot=F)
tz <-hist(z1,breaks=br,plot=F)
ylim <- max(tx\$count,ty\$count,tz\$count)+500
xvalue=tx\$mid
plot(xvalue,tx\$count,ylim=c(0,ylim),xlim=c(0,xlim),main=title,xlab=xlab,ylab=ylab,type="b",pch=1,col=2)
lines(xvalue,ty\$count,ylim=c(0,ylim),xlim=c(0,xlim),type="b",pch=2,col=3)
lines(xvalue,tz\$count,ylim=c(0,ylim),xlim=c(0,xlim),type="b",pch=3,col=4)
legend("topright",c("share","unique","noncluster"),lty=1,pch=1:3,lwd=1,col=c(2,3,4))
}

pdf("$title.pdf")

plotlengthline(200,x1[,2],y1[,2],z1[,2],5000,"cds_length","Length (bp)","Count")
plotlengthline(200,x2[,2],y2[,2],z2[,2],1000,"exon_length","Length (bp)","Count")
plotnumberline(1,x3[,2],y3[,2],z3[,2],10,"exon_number","Exon Number","Count")
plotlengthline(200,x4[,2],y4[,2],z4[,2],1000,"intron_length","Length (bp)","Count")
plotlengthline(200,x5[,2],y5[,2],z5[,2],5000,"mRNA_length","Length (bp)","Count")
plotlengthline(1000,x6[,2],y6[,2],z6[,2],50000,"intergenic_length","Length (bp)","Count")
dev.off()
END.
close OUT;
system ("cat $title.r | R --vanilla --slave");
}
 
