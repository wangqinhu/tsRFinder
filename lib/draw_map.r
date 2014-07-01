f<-list.files("tmp")
for (i in 1:length(f)) {
  n<-unlist(strsplit(f[i], "[.]"))[1]
  pdf(paste("img/", n, ".pdf", sep=""))
  data<-read.table(paste("tmp/", f[i], sep=""), header=TRUE)
  mex<-max(data$plus,data$minus)
  plot(data$position, data$plus, type="h", col="green", xlab="Offset", ylab="Counts", ylim=range(-mex:mex))
  abline(h=0)
  points(data$position, -data$minus, type="h", col="black")
  dev.off()
}
rm(list=ls())
