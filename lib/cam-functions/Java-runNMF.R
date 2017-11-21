
# This software is under license BSD 3.
# ------------------ Perform NMF Algorithm ---------------------------------------

packageExist <- require("NMF")
if (!packageExist) {
  install.packages('NMF', repos=c('http://web.cbio.uct.ac.za/~renaud/CRAN', getOption('repos')))
  library("NMF")
}

cat("Performing NMF to decompose data ... \n")
      
NMFresult <- nmf(X_mask,K)
A_est <- NMFresult@fit@H
A_est <- t(A_est)
S_est <- NMFresult@fit@W

cat("Done\n")

# pass results to Java environment
final.result <- list(Aest=A_est)

# save data
source("Java-saveData.R")
saveData(A_est, S_est)


detach()