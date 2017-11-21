
# This software is under license BSD 3.

# ------------------ Perform fastICA Algorithm ---------------------------------------
# load the library "fastICA" for the function "fastICA"
packageExist <- require("fastICA")
if (!packageExist) {
  install.packages("fastICA")
  library("fastICA")
}

cat("Performing fastICA to decompose data ... \n")
                 
ICAresult <- fastICA(X_mask,K)
A_est <- t(ICAresult$A)
S_est <- ICAresult$S

cat("Done\n")

# pass results to Java environment
final.result <- list(Aest=A_est)

# save data
source("Java-saveData.R")
saveData(A_est, S_est)


detach()
