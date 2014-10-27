pdf(file="distribution.pdf", 6, 9)
# read length file
rls<-read.table("srna.len")
rlt<-read.table("trna.len")
rld<-rlt
rld[1:length(rls$V1),3]<-rls$V2 - rlt$V2
# plot
layout(c(1:3))
par(mar=(c(4.5,5,2,1)))
barplot(rls$V2, col=2, main="Length distribution of sRNA reads", xlab="length (nt)", ylab="Frequency", names.arg=rls$V1)
par(mar=(c(4.5,5,2,1)))
barplot(rlt$V2, col=3, main="Length distribution of tRNA reads", xlab="length (nt)", ylab="Frequency", names.arg=rlt$V1)
par(mar=(c(4.5,5,2,1)))
barplot(t(rld), col=c("red","green"), main="Length distribution of tRNA and non-tRNA reads", xlab="length (nt)", ylab="Frequency", names.arg=rlt$V1)
dev.off()
