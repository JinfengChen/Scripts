
pvalue <- function(x1, y1, x2, y2, top, p){
     star <- '*'
     if (p > 0.1) {star <- ''}
     if (p < 0.001){ star <- '**'}
     if (p < 0.0001){ star <- '***'}
     segments(x1,y1,x1,top)
     segments(x1,top,x2,top)
     segments(x2,top,x2,y2)
     #segments(x1-0.2,y1,x1+0.2,y1)
     #segments(x2-0.2,y2,x2+0.2,y2)
     xt <- min(x1,x2)+abs(x2-x1)/2
     yt <- top + 100
     text(xt,yt,star, cex=1.5)
} 

pdf("mPing_intron_distance_boxplot.pdf")

som <- read.table("Somatic.intron.distance")
str <- read.table("Strains.intron.distance")
ril <- read.table("RIL.intron.distance")
sim <- read.table("Simulation.intron.distance")

y <- wilcox.test(str[,1],som[,1],alternative="greater", correct=FALSE)
x <- wilcox.test(str[,1],sim[,1],alternative="greater", correct=FALSE)

boxplot(som[,1], ril[,1], str[,1], sim[,1], ylim=c(0,4000), outline=FALSE, col=c("aquamarine3", "steelblue2" ,"sandybrown","dim gray"), ylab="Length of Intron (bp)",xaxt='n',frame.plot=FALSE)
#legend('topright', bty='n', border='NA', lty=c(0,0),cex=1 ,c("Somatic", "Strains", "RIL", "Simulation"),fill=c("aquamarine3", "steelblue2" ,"sandybrown", "dim gray"))
text(1:4,rep(-500,7), cex=1, offset=2,labels=c("Somatic", "RIL", "Strains", "Simulation"),srt=55,xpd=TRUE)
axis(1,c(0.5,4.5), labels=c("",""),line=0)
pvalue(3.1,3600,4,2700,3700,x$p.value)
pvalue(2.9,3600,1,3500,3700,y$p.value)

dev.off()

