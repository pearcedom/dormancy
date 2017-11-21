# This software is under license BSD 3.
# ----------------------------- Perform CAM Algorithm ---------------------------------------
source("CAM.R")
algorithm  <- 0                         # indicate CAM-CM
CAM.result <- CAM(t(X_mask),K,denoise,vis) # all parameters will be preloaded by Java app.
A_est <- CAM.result[[1]]
S_est <- CAM.result[[2]]
# ---------------------------------- End of CAM ---------------------------------------------


# Decide which curve corresponds to the input function
peak   <- rep(0,K)
peakIX <- rep(0,K)
for (i in 1:K) {
  peak[i]   <- max(A_est[,i])
  peakIX[i] <- which.max(A_est[,i])
}

index.input <- which.max(peak) # the input function has largest peak
index.other <- (1:K)[-index.input]
tmp <- sort(peakIX[index.other], index.return=TRUE)
index.other <- index.other[tmp$ix]

TC_p     <- A_est[,index.input,drop=FALSE]
TC_other <- A_est[,index.other,drop=FALSE]
TCs <- cbind(TC_other, TC_p)


# ----------------------------- Perform CM Algorithm ----------------------------------------
# random initialization
initk <- matrix(0,K,2)
for (i in 1:K) {
  initk[i,1] <- 0.1 * runif(1)
  initk[i,2] <- initk[i,1] + runif(1)
}

source("CM.R")
CM.result <- CM(t(X_mask),TCs,initk,del_t)
eKtrans         <- CM.result[[1]]
eKep            <- CM.result[[2]]
pixelwise_Ktran <- CM.result[[3]]

kinetic   <- matrix(0, 2, K-1)
for (i in 1:(K-1)) {
  kinetic[1,i] <- eKtrans[i]
  kinetic[2,i] <- eKep[i]
}
dimnames(kinetic)[[1]] <- paste(c("K_in","K_out"), "(/min)")
dimnames(kinetic)[[2]] <- paste("Compartment", 1:(K-1), sep="")
cat("\nThe CAM-CM estimated kinetic parameters:\n")
print(kinetic)

final.result <- list(Aest=TCs, cmResults=kinetic)
# ---------------------------------- End of CM ----------------------------------------------


# save data
source("Java-saveData.R")
saveData(TCs, t(S_est))


detach()

