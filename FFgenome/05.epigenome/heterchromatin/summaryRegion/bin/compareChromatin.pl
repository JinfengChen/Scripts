#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"inf1:s","inf2:s","project:s","help");


my $help=<<USAGE;
Compare information of gene length/exon, intron length/intergenic length/gene expression/methylation/flanking TE of gene for euchromatin and heterchromatin for two species. 
perl $0 --inf1 OBa_manual_chromatin --inf2 rice_chromatin --project OBa2rice > log 2> log2 &
--inf1: chromatin information for species1
--inf2: chromatin information for species2
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

`mkdir $opt{project}` unless (-f $opt{project});

####species 1 
my $eu1="$opt{inf1}/euchromatin";
my $he1="$opt{inf1}/heterochromatin";
my $eu2="$opt{inf2}/euchromatin";
my $he2="$opt{inf2}/heterochromatin";
####gene structure
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"mRNA_length","Length (bp)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"cds_length","Length (bp)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"exon_length","Length (bp)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"intron_length","Length (bp)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"intron1_length","Length (bp)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"intron2_length","Length (bp)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"exon_number","Number");
####gene expression
& compareStructureDrawBar4Expression($eu1,$he1,$eu2,$he2,"expression","Normalized expression level");
#### TE
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"TEdensity","Density (0-1)");
#### methylation
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"Upstream","Methylation Level (CG)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"Body","Methylation Level (CG)");
& compareStructureDrawBar($eu1,$he1,$eu2,$he2,"Downstream","Methylation Level (CG)");

#######################
sub compareStructureDrawLine
{
my ($eu,$he,$title)=@_;
open OUT, ">$opt{project}/$title.r" or die "$!";
print OUT <<"END.";
read.table("$eu/$title.txt") -> x
read.table("$he/$title.txt") -> y
pdf("$opt{project}/$title.pdf")
plot(density(x[,2]),col=2);
lines(density(y[,2]),col=3);
dev.off()
END.
close OUT;
system ("cat $opt{project}/$title.r | R --vanilla --slave");
}

sub compareStructureDrawBar
{
my ($eu1,$he1,$eu2,$he2,$title,$unit)=@_;
open OUT, ">$opt{project}/$title.r" or die "$!";
print OUT <<"END.";
read.table("$eu1/$title.txt") -> eu1
read.table("$he1/$title.txt") -> he1
read.table("$eu2/$title.txt") -> eu2
read.table("$he2/$title.txt") -> he2
x <- eu1[,2]
y <- he1[,2]
a <- eu2[,2]
b <- he2[,2]
name <-c("euchromatin","heterochromatin")
##
xmean <- mean(x)
xse <- 1.96*sd(x)/sqrt(length(x))
xcilow <- xmean-xse
xcihigh <- xmean+xse

##
ymean <- mean(y)
yse <- 1.96*sd(y)/sqrt(length(y))
ycilow <- ymean-yse
ycihigh <- ymean+yse

##
amean <- mean(a)
ase <- 1.96*sd(a)/sqrt(length(a))
acilow <- amean-ase
acihigh <- amean+ase 

##
bmean <- mean(b)
bse <- 1.96*sd(b)/sqrt(length(b))
bcilow <- bmean-bse
bcihigh <- bmean+bse

##
spe1mean <- c(xmean,ymean)
spe2mean <- c(amean,bmean)
spe1cilow <- c(xcilow,ycilow)
spe2cilow <- c(acilow,bcilow)
spe1cihigh <- c(xcihigh,ycihigh)
spe2cihigh <- c(acihigh,bcihigh)


pdf("$opt{project}/$title.pdf",width=4,height=4)
mean <- rbind(spe1mean,spe2mean)
cilow <- rbind(spe1cilow,spe2cilow)
cihigh <- rbind(spe1cihigh,spe2cihigh)
topvalue=max(mean)+0.4*max(mean)
legendy=max(mean)+0.2*max(mean)
library("gplots")
barplot2(mean,ylim=c(0,topvalue),main="$title",ylab="$unit",beside = TRUE,plot.ci=TRUE,col = c("gray38", "gray57"),names.arg=name,ci.l=cilow,ci.u=cihigh)
savefont <- par(font=3)
legend(1,topvalue,bty="n",c("O. brachyantha","O. sativa"),lty=1,lwd=12,col=c("gray38", "gray57"))
abline(h=0)
#text(barx,par("usr")[3],name,srt=0,xpd=TRUE)
####lineplot
plot(density(x,n=50),cex=0.4,type="b",col=2,pch=1,xlab="$unit",ylab="Density",main="$title");
lines(density(y,n=50),cex=0.4,type="b",col=2,pch=2);
lines(density(a,n=50),cex=0.4,type="b",col=3,pch=3);
lines(density(b,n=50),cex=0.4,type="b",col=3,pch=4);
savefont <- par(font=3)
legend("topright",cex=0.5,c("O. brachyantha-EU","O. brachyantha-HE","O. sativa-EU","O. sativa-HE"),pch=c(1:4),lty=c(1,1,1,1),col=c(2,2,3,3))
dev.off()
END.
close OUT;
system ("cat $opt{project}/$title.r | R --vanilla --slave");
}

sub compareStructureDrawBar4Expression
{
my ($eu1,$he1,$eu2,$he2,$title,$unit)=@_;
open OUT, ">$opt{project}/$title.r" or die "$!";
print OUT <<"END.";
read.table("$eu1/$title.txt") -> eu1
read.table("$he1/$title.txt") -> he1
read.table("$eu2/$title.txt") -> eu2
read.table("$he2/$title.txt") -> he2
x <- eu1[,2]
y <- he1[,2]
a <- eu2[,2]
b <- he2[,2]
name <-c("euchromatin","heterochromatin")
##
xmean <- mean(x)
xse <- 1.96*sd(x)/sqrt(length(x))
xcilow <- xmean-xse
xcihigh <- xmean+xse

##
ymean <- mean(y)
yse <- 1.96*sd(y)/sqrt(length(y))
ycilow <- ymean-yse
ycihigh <- ymean+yse

##
amean <- mean(a)
ase <- 1.96*sd(a)/sqrt(length(a))
acilow <- amean-ase
acihigh <- amean+ase 

##
bmean <- mean(b)
bse <- 1.96*sd(b)/sqrt(length(b))
bcilow <- bmean-bse
bcihigh <- bmean+bse

##
spe1mean <- c(xmean,ymean)
spe2mean <- c(amean,bmean)
spe1cilow <- c(xcilow,ycilow)
spe2cilow <- c(acilow,bcilow)
spe1cihigh <- c(xcihigh,ycihigh)
spe2cihigh <- c(acihigh,bcihigh)


pdf("$opt{project}/$title.pdf",width=4,height=4)
mean <- rbind(spe1mean,spe2mean)
cilow <- rbind(spe1cilow,spe2cilow)
cihigh <- rbind(spe1cihigh,spe2cihigh)
####barplot
library("gplots")
barplot2(mean,ylim=c(-0.12,0.12),main="$title",ylab="$unit",beside = TRUE,plot.ci=TRUE,col = c("gray38", "gray57"),names.arg=name,ci.l=cilow,ci.u=cihigh)
savefont <- par(font=3)
legend(1,0.12,c("O. brachyantha","O. sativa"),bty="n",lty=1,lwd=8,col=c("gray38", "gray57"))
abline(h=0)
#text(barx,par("usr")[3],name,srt=0,xpd=TRUE)
####lineplot
plot(density(x,n=50),cex=0.4,type="b",col=2,pch=1,xlab="$unit",ylab="Density",main="$title");
lines(density(y,n=50),cex=0.4,type="b",col=2,pch=2);
lines(density(a,n=50),cex=0.4,type="b",col=3,pch=3);
lines(density(b,n=50),cex=0.4,type="b",col=3,pch=4);
savefont <- par(font=3)
legend("topleft",cex=0.5,c("O. brachyantha-EU","O. brachyantha-HE","O. sativa-EU","O. sativa-HE"),pch=c(1:4),lty=c(1,1,1,1),col=c(2,2,3,3))
dev.off()
END.
close OUT;
system ("cat $opt{project}/$title.r | R --vanilla --slave");
}


=pod
sub compareStructure11
{
open OUT, ">$opt{project}.r" or die "$!";
print OUT <<"END.";
read.table("$opt{project}.4r") -> x
pdf("$opt{project}.pdf")
library("gplots")
barplot2(x[,2],plot.ci=TRUE,ci.l=x[,3],ci.u=x[,4],names.arg=x[,1],ylab="Normalized expression difference")
mean <- rbind(x[,6],x[,10])
cilow <- rbind(x[,7],x[,11])
cihigh <- rbind(x[,8],x[,12])
barplot2(mean,ylim=c(0,1.5),beside = TRUE,plot.ci=TRUE,col = c("gray38", "gray57"),names.arg=x[,1],ci.l=cilow,ci.u=cihigh,y
savefont <- par(font=3)
legend(1,1.4,c("O. brachyantha","O. sativa"),lty=1,lwd=12,col=c("gray38", "gray57"))
abline(h=0)
dev.off()
END.
close OUT;
q
system ("cat $opt{project}.r | R --vanilla --slave");
}
=cut


sub ci
{
my ($num)=@_;
my $loop=0;
my $total;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}
my $mean=$total/$number;
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return ($mean,$se,$number);
}

