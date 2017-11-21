# This software is under license BSD 3.
# Normalize TCs
TCs.norm <- t(t(TCs)/colSums(TCs))


# Curve fitting by estimated kinetic parameters (CM results)
L   <- nrow(TCs)
xx  <- c(4.5, 1.2, 1.5, 0.005, -1.1)
eCp <- rep(0,L)
for (t in 1:L) {
  eCp[t] <- xx[1]*exp(-xx[2]*t) + xx[3]*exp(-xx[4]*t) + xx[5]
}
eCp <- c(TCs[L,K], eCp[2:L])

eC  <- matrix(0,L,K)
for (i in 1:(K-1)) {
  ki <- eKtrans[i]; ko <- eKep[i]
  H  <- matrix(0,L,L)
  for (t in 1:L) {
    H[t,1]   <- ki * exp(-(t-1) * del_t * ko)
    if (t < L) {
      for (tt in 1:(L-t))
        H[t+tt,tt+1] <- H[t,1]
    }
  }
  eC[,i] <- H %*% eCp
}
eC[,K] <- eCp
eC.norm <- t(t(eC)/colSums(eC))


# ---------------- Draw Figures of Estimated Contrast/Tracer Concentration) -----------------
x <- (1:nrow(TCs.norm)) * del_t
xmax  <- max(x)
ymax  <- ceiling(max(c(TCs.norm,eC.norm)) / 0.02) * 0.02

plot.new()
par(mar=c(5,5,4,2), tck=-0.02, tcl=NA, las=1, font.lab=2, font.axis=2)
par(usr=c(0,xmax,0,ymax))

axis(1, at=0:xmax,             col="blue", col.axis="blue")
axis(2, at=seq(0,ymax,0.02),   col="blue", col.axis="blue")
lines(x, TCs.norm[,K], type="p", lwd=2, col="blue")
lines(x, TCs.norm[,1], type="p", lwd=2, col="red")
lines(x, TCs.norm[,2], type="p", lwd=2, col="green")
lines(x, eC.norm[,K], lty="solid",   lwd=2, col="blue")
lines(x, eC.norm[,1], lty="dashed",  lwd=2, col="red")
lines(x, eC.norm[,2], lty="dotdash", lwd=2, col="green")
title(xlab="Time (minutes)", ylab="Normalized TC", col.lab="blue", 
      main="Estimated Contrast/Tracer Concentration")
legend("topright", c("input","flow1","flow2"),
       lty=c("solid","dashed","dotdash"), lwd=2, col=c("blue","red","green"))
# -------------------------------- End of Drawing -------------------------------------------