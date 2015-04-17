pdf("unalign.anno.pdf")
read.table("HEG4_RAW.unalign.summary",skip=1) -> ref
lbls <- paste(ref[,1],floor(ref[,3]*100))
lbls <- paste(lbls,"%",sep="")
pie(ref[,2],labels=lbls,radius=0.8,col=rainbow(length(lbls))) -> x
#legend(1,0.5,cex=0.6,ref[,1],bty="n")

read.table("HEG4_luluv2.unalign.summary",skip=1) -> ref
lbls <- paste(ref[,1],floor(ref[,3]*100))
lbls <- paste(lbls,"%",sep="")
pie(ref[,2],labels=lbls,radius=0.8,col=rainbow(length(lbls))) -> x

read.table("HEG4_base_uncover.summary",skip=1) -> ref
lbls <- paste(ref[,1],floor(ref[,3]*100))
lbls <- paste(lbls,"%",sep="")
pie(ref[,2],labels=lbls,radius=0.8,col=rainbow(length(lbls))) -> x



dev.off()

