# Read data
tsr5<-read.table("tscs.5")
tsr3<-read.table("tscs.3")

# min-max-range
xmin<-min(tsr5$V1,tsr3$V1)
xmax<-max(tsr5$V1,tsr3$V1)
ymin<-min(tsr3$V2)
ymax<-max(tsr5$V2)
xr<-seq(xmin-1,xmax,1)
xleft<-xmin-1
xright<-xmax+1
ybottom<-ymin-2
ytop<-ymax+2

pdf("cleavage_profile.pdf", 8, 4)
par(mar=c(4,1,1,0))
# plot area
plot(c(xleft, xright), c(ybottom, ytop), type = "n",
     xlab = "tRNA base", ylab = "",
     main="Distribution of tRNA clevage sites",
     axes=F
     )
axis(side=1, labels=T,  at = seq(xmin-1,xmax+1))
mtext(side = 2, text = "Cleavage frequency", line = 0)

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
i<-1:ymax
rect(xleft, 1+i, xleft+0.4, 1+i+1, col=topo.colors(max(tsr5$V2)), border = "white")

# -y
j<-ymin:-1
rect(xleft, j-1, xleft+0.4, j, col=rev(cm.colors(max(-tsr3$V2))), border = "white")
dev.off()