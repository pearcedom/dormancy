# This software is under license BSD 3.

ltype <- c("solid","dashed","dotdash","dotted","longdash","twodash")
color <- c("blue","red","green","brown","pink","violet")


# Normalize TCs
TCs.norm <- t(t(TCs)/colSums(TCs))

x <- (1:nrow(TCs.norm)) * del_t
xmax  <- max(x)
if (algorithm == 0) {
  ymax  <- ceiling(max(TCs.norm) / 0.02) * 0.02
} else {
  ymax  <- ceiling(max(TCs.norm) / 0.2) * 0.2
}


# ---------------- Draw Figures of Estimated Contrast/Tracer Concentration) -----------------
plot.new()
par(mar=c(5,5,4,2), tck=-0.02, tcl=NA, las=1, font.lab=2, font.axis=2)
par(usr=c(0,xmax,0,ymax))

if (algorithm == 0) {
  axis(1, at=0:xmax,             col="blue", col.axis="blue")
  axis(2, at=seq(0,ymax,0.02),   col="blue", col.axis="blue")
  
  lines(x, TCs.norm[,K], lty=ltype[1], lwd=2, col=color[1])
  for (i in 1:(K-1))
    lines(x, TCs.norm[,i], lty=ltype[i%%6+1],  lwd=2, col=color[i%%6+1])
  
  title(xlab="Time (minutes)", ylab="Normalized TC", col.lab="blue", 
        main="Estimated Contrast/Tracer Concentration")
  legend("topright", c("input",paste("flow",1:(K-1),sep="")),
        lty=ltype[(1:K-1)%%6+1], lwd=2, col=color[(1:K-1)%%6+1])
} else {
  axis(1, at=seq(0,xmax,del_t), col="blue", col.axis="blue")
  axis(2, at=seq(0,ymax,0.2),   col="blue", col.axis="blue")
  
  for (i in 1:K)
    lines(x, TCs.norm[,i], lty=ltype[(i-1)%%6+1], lwd=2, col=color[(i-1)%%6+1])
  
  title(xlab="Time (minutes)", ylab="Normalized TC", col.lab="blue", 
        main="Estimated Contrast/Tracer Concentration")
  legend("topright", paste("col",1:K,sep=""),
        lty=ltype[(1:K-1)%%6+1], lwd=2, col=color[(1:K-1)%%6+1])
}
# -------------------------------- End of Drawing -------------------------------------------