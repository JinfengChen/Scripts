require(Hmisc)
pdf("indel_codon_len.length_indel_length.pdf")

par(mar=c(6,4,4,2))
x <- read.table("indel_codon_len.length.distr",skip=1)
data <- rbind(x[,3]/sum(x[,3]),x[,4]/sum(x[,4]))
xx <- barplot(data,beside=TRUE,ylab="Proportion",border=FALSE,ylim=c(0,0.7),col=c("Orange","blue"))

x1 <- read.table("indel.length.distr",skip=1)
data1 <- rbind(x1[,3]/sum(x1[,3]),x1[,4]/sum(x1[,4]))
tmp2=subplot(xx1 <- barplot(data1,beside=TRUE,border=FALSE,ylim=c(0,0.7),cex=0.5,cex.axis=0.7,cex.lab=0.7,tck=-0.02,col=c("Orange","blue"),mgp=c(0.8,0.8,0)),30,0.5,c(3,3))
op <- par(no.readonly=TRUE)
par(tmp2)
axis(1,c(0.5,max(xx1)+0.5),line=0,labels=c("",""))
text(xx1[1,]+0.2,rep(-0.02,6),offset=2,cex=0.7,tck=-0.02,labels=x1[,2],srt=55,xpd=TRUE)
par(op)


axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
text(xx[1,]+0.2,rep(-0.02,6),offset=2,labels=x[,2],srt=55,xpd=TRUE)
text(xx[1,],x[,3]/sum(x[,4])+0.03,offset=2,labels=x[,3],srt=55,xpd=TRUE)
text(xx[1,]+1.2,x[,4]/sum(x[,4])+0.03,offset=2,labels=x[,4],srt=55,xpd=TRUE)
legend("topright",c("DEL","INS"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("orange","blue"))
mtext("INDEL Length (bp)",side=1, at=23,line=3)
dev.off()

