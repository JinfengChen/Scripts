error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pdf("ping_number.pdf")
data =read.table("RIL275_RelocaTEi.CombinedGFF.characterized.ping_number.summary", header = T)
expr = data[,2]
std = data[,3]
barx <- barplot(expr, col=c("black"), ylim=c(0,120), border=F, axis.lty=1, xlab='', ylab='')
error.bar(barx, expr, std)
axis(1,c(0.1, max(barx)+0.6),line=0,labels=c("",""))
text(barx, rep(-5, 6),offset=2,labels=data[,1],srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 3.6
ypos <- 44
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Number", side=1,font=1, at=xpos+1.1,line=3, cex=1, col="black")
mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
mtext("mPing", side=2,font=3, at=ypos+11,line=3, cex=1, col="black")
mtext("number", side=2,font=1, at=ypos+25,line=3, cex=1, col="black")
dev.off()


