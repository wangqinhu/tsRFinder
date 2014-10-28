pdf(file="distribution.pdf", 6, 9)
# read length file
rls<-read.table("srna.len")
rlt<-read.table("trna.len")
len<-length(rls$V1)
rld=rep(0, len*3)
dim(rld)<-c(len,3)
rld<-rlt[,1:2]
rld[1:len,3]<-rls$V2 - rlt$V2
# plot
layout(c(1:3))
par(mar=(c(4.5,5,2,1)))
barplot(t(rls[,3:6]), col=c("green","red","blue","black"), legend = c("A","T","C","G"), main="Length distribution of sRNA reads", xlab="length (nt)", ylab="Frequency", names.arg=rls$V1)
par(mar=(c(4.5,5,2,1)))
barplot(t(rlt[,3:6]), col=c("green","red","blue","black"), legend = c("A","T","C","G"), main="Length distribution of tRNA reads", xlab="length (nt)", ylab="Frequency", names.arg=rlt$V1)
par(mar=(c(4.5,5,2,1)))
barplot(t(rld[,1:3]), col=c("red","green"), legend = c("sRNA","tRNA"), main="Length distribution of tRNA and non-tRNA reads", xlab="length (nt)", ylab="Frequency", names.arg=rlt$V1)
dev.off()
