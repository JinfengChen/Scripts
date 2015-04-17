pdf("HEG4_luluv2.CDS.uncover.overlap.missrate.4r.pdf")
par(mai=c(1.7,1,1,0.5))
read.table("HEG4_luluv2.CDS.uncover.overlap.missrate.4r") -> x

barplot(x[,3],axes=FALSE,border=F,ylim=c(0,2000),col=c("cornflowerblue")) -> xx
for (i in 1:length(xx)) { # adj can not take vector, so we use loops to add text
  text(xx[i],-50,labels=x[i,2],cex=1,srt=65,adj=c(1,1),xpd=TRUE)
}
axis(1,c(0,max(xx)+0.7),labels=c("",""))
axis(2,seq(0,2000,by=2000/4),cex=1.2)
text(min(xx)-2.5,2000/2,labels="Number of Gene",srt=90,cex=1.1,xpd=TRUE)
text(max(xx)/2+0.3,-2000/4,labels="Percent of Gap (%)",cex=1.1,xpd=TRUE)

dev.off()

