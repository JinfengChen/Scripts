pdf("lengthclass.Qry.pdf")
par(mar=c(6,4,4,10))
par(xpd=TRUE)

barlegend <- function(x,y,height,length,name,color){
    rect(x,y-0.04*height,x+0.05*length,y,col=color,border=FALSE)
    text(x+0.15*length,y-0.02*height,labels=name)
}


x <- read.table("Qry.indel.anno.lengthclass")
data <- rbind(x[,3]/x[,2],x[,4]/x[,2],x[,5]/x[,2],x[,6]/x[,2],x[,7]/x[,2],x[,8]/x[,2])
xx <- barplot(data,ylab="Proportion",border=FALSE,space=1,ylim=c(0,1),col=c("red","orange","blue","green","gray","black"))
axis(1,c(0.9,max(xx)+0.6),line=0,labels=c("",""))
text(xx,rep(-0.07,6),offset=2,labels=x[,1],srt=55,xpd=TRUE)
legend(8.5,0.8,bty="n",lty=c(0,0),border="NA",cex=1.2,c("CDS","intron","LTR","DNA TE","other TE","unknown"),fill=c("red","orange","blue","green","gray","black"))

#barlegend(10,0.7,1,10,"CDS","red")

mtext("Indel Length (bp)",side=1, at=4.8,line=4)
dev.off()


