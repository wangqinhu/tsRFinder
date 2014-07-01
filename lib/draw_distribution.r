pdf(file="distribution.pdf", 6, 6)
# read length file
rls<-read.table("srna.len")
rlt<-read.table("trna.len")
# plot
layout(c(1,2))
par(mar=(c(4.5,5,2,1)))
barplot(rls$V2, col=2, main="sRNA read", xlab="length (nt)", ylab="Frequency", border="white", names.arg=rls$V1)
par(mar=(c(4.5,5,2,1)))
barplot(rlt$V2, col=3, main="tRNA read", xlab="length (nt)", ylab="Frequency", border="white", names.arg=rlt$V1)
dev.off()
