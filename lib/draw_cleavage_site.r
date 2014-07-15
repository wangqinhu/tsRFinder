# Read data
tsr5<-read.table("tscs.5")
tsr3<-read.table("tscs.3")

# min-max-range
xmin<-min(tsr5$V1,tsr3$V1)
xmax<-max(tsr5$V1,tsr3$V1)
yu<-max(tsr5$V2)
yd<-min(tsr3$V2)
xr<-seq(xmin-1,xmax,1)

pdf("cleavage_profile.pdf", 8, 4)
# plot area
plot(c(xmin, xmax), c(yd-5, yu+5), type = "n", xlab = "tRNA base", ylab = "Cleavage frequency", main="Distribution of tRNA clevage sites", axes=F)

# tRNA
rect(xr, 0, xr+0.8, 1, col = 1, lwd = 0.5)

# anti codon
ac<-c(0,1,2)
rect(ac, 0, ac+0.8, 1, col = 2, lwd = 0.5)

# cleavage
# 5' tsR
rect(tsr5$V1 - 0.4,  2, tsr5$V1 + 0.4, tsr5$V2 + 2,
     border = "white", 
     col= topo.colors( max( tsr5$V2) )[ tsr5$V2 ]
)
# 3' tsR
rect(tsr3$V1 - 0.4, -1, tsr3$V1 + 0.4, tsr3$V2 - 1,
     border = "white", 
     col= cm.colors( max( -tsr3$V2) )[ -tsr3$V2 ]
)

# define axis
# +y
i<-1:max(tsr5$V2)
rect(xmin-1.4, 1+i , xmin-1.1, 1+i+1, col=topo.colors(max(tsr5$V2)), border = "white")

# -y
j<-min(tsr3$V2):-1
rect(xmin-1.4, j-1,    xmin-1.1, j, col=rev(cm.colors(max(-tsr3$V2))), border = "white")
dev.off()
