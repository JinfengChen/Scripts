error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pdf("Qpcr.pdf")
data =read.table("20140323.sum.txt", header = T)
expr = rbind(data[,2], data[,4])
std = rbind(data[,3], data[,5])
barx <- barplot(expr, beside=TRUE, col=c("blue", "orange"), ylim=c(0,2), border=F, axis.lty=1, xlab="Gene", ylab="Relative Expression Level")
error.bar(barx, expr, std)
axis(1,c(0.9,max(barx)+0.6),line=0,labels=c("",""))
text(barx[1,]+0.5,rep(-0.07,6),offset=2,labels=data$Gene,srt=0,xpd=TRUE)
legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
dev.off()


