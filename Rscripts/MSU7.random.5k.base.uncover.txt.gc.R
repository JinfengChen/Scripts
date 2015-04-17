pdf("MSU7.random.5k.base.uncover.txt.gc.pdf")
read.table("MSU7.random.5k.base.uncover.txt.gc") -> x

require(Hmisc)
plot(x[,2],x[,3],xlim=c(0,1000),xlab="Length (bp)",ylab="GC (%)")
subplot(plot(x[,2],x[,3],xlab="Length (bp)",ylab="GC (%)",cex=0.5,cex.axis=0.5,cex.lab=0.5,tck=-0.02,mgp=c(0.8,0.3,0)),800,0.15)
#axis(1,c(0,max(xx)+0.7),labels=c("",""))
#axis(2,seq(0,80,by=20),cex=1.2)
#text(min(xx)-2,40,labels="Percent of Genome (%)",srt=90,cex=1.1,xpd=TRUE)
#text(max(xx)/2+0.3,-20,labels="Sequence Depth",cex=1.1,xpd=TRUE)
dev.off()

