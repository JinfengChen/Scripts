pdf("Large_effect.Pfam.summary.4test.summary.Q0.05.pdf", height=6, width=8)

par(mar=c(6,24,4,2))
x <- read.table("Large_effect.Pfam.summary.4test.summary.Q0.05",sep="\t")
data <- x[,2]
xx <- barplot(data,space=0.5,horiz=TRUE,xlab="# of genes with large-effect mutations",border=FALSE,col=c("cornflowerblue"))
axis(2,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
text(rep(-3,7),xx[,1],adj=1, offset=2,labels=x[,9],xpd=TRUE)
#text(xx[1,],x[,4]/sum(x[,4])+0.03,offset=2,labels=x[,4],srt=55,xpd=TRUE)
#text(xx[1,]+1.2,x[,5]/sum(x[,5])+0.03,offset=2,labels=x[,5],srt=55,xpd=TRUE)
#legend(15,0.7,bty="n",lty=c(0,0),cex=1.2,c("Deletion","Insertion"),fill=c("Orange","blue"))
#barlegend(11,0.47,0.6,16,"Deletion","Orange")
#barlegend(11,0.44,0.6,16,"Insertion","Blue")
#mtext("Indel Length (bp)",side=1, at=7,line=4)

dev.off()

