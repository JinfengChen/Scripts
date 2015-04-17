#pdf("coverage.curve.pdf")
#read.table("MSU7.100bpwin.coverage") -> x
#pdf("coverage.random.curve.pdf")
#read.table("MSU7.random.100bpwin.coverage") -> x
#pdf("coverage.random.5k.curve.pdf")
#read.table("MSU7.random.5kwin.coverage") -> x
pdf("coverage.5k.curve.pdf")
read.table("MSU7.5kwin.coverage") -> x
a <- c(seq(0,400,5),50000);
hist(x[,4],breaks=a,plot = FALSE) -> xh
atx <- seq(0,400,50)
aty <- seq(0,0.016,0.004)
plot(xh$mids,xh$density,xlim=c(0,400),ylim=c(0,0.016),type="p",pch=18,col="blue",xlab="Coverage",ylab="Frequency",axes=FALSE);
axis(1,at=atx,labels=atx)
axis(2,at=aty,labels=aty)
dev.off()

