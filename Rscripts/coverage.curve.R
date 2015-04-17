pdf("coverage.curve.pdf")
read.table("MSU7.100bpwin.coverage") -> x
hist(x[,4],plot = FALSE) -> xh
atxvalue <- seq (0,max(xh$breaks),by=20)
atyvalue <- seq (0,max(xh$density)*1000,by=10)
atx <- c(atxvalue,max(atxvalue)+20)
aty <- c(atyvalue,max(atyvalue)+10)
plot(xh$mids,1000*xh$density,type="p",pch=18,col="blue",xlab="Coverage",ylab="Frequency",xlim=c(0,max(atxvalue)+20),ylim=c(0,max(atyvalue)+10),axes=FALSE)
axis(1,at=atx,labels=atx)
axis(2,at=aty,labels=aty)
dev.off()

