pdf("length.pdf")

barlegend <- function(x,y,height,length,name,color){
    rect(x,y-0.04*height,x+0.05*length,y,col=color,border=FALSE)
    text(x+0.15*length,y-0.02*height,labels=name)
}


par(mar=c(6,4,4,2))
x <- read.table("length.distr",skip=1)
data <- rbind(x[,4]/sum(x[,4]),x[,5]/sum(x[,5]))
xx <- barplot(data,beside=TRUE,ylab="Proportion",border=FALSE,ylim=c(0,0.5),col=c("Orange","blue"))
axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
text(xx[1,]+0.2,rep(-0.04,6),offset=2,labels=x[,2],srt=55,xpd=TRUE)
text(xx[1,],x[,4]/sum(x[,4])+0.03,offset=2,labels=x[,4],srt=55,xpd=TRUE)
text(xx[1,]+1.2,x[,5]/sum(x[,5])+0.03,offset=2,labels=x[,5],srt=55,xpd=TRUE)
#legend(15,0.7,bty="n",lty=c(0,0),cex=1.2,c("Deletion","Insertion"),fill=c("Orange","blue"))
barlegend(11,0.47,0.6,16,"Deletion","Orange")
barlegend(11,0.44,0.6,16,"Insertion","Blue")
mtext("Indel Length (bp)",side=1, at=7,line=4)



dev.off()

