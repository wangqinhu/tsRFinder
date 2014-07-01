x<-read.table("infile.txt", header=TRUE)
a<-x$plus
b<-floor(9*(a-min(a))/(max(a)-min(a)))
write.table(b, "BDI.txt", quote=FALSE, row.names=FALSE,col.names=FALSE)
