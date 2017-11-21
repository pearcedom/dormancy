
# This software is under license BSD 3.

# ------------------ Perform Factor Analysis Algorithm ---------------------------------------
packageExist <- require("MASS")
if (!packageExist) {
  install.packages("MASS")
  library("MASS")
}

# load the library "stats" for the function "factanal"
packageExist <- require("stats")
if (!packageExist) {
  install.packages("stats")
  library("stats")
}

cat("Performing Factor Analysis to decompose data ... \n")

FAresult <- factanal(X_mask,K,scores="regression")
A_est <- matrix((FAresult$loadings),ncol=K)
S_est <- FAresult$scores

cat("Done\n")

# pass results to Java environment

final.result <- list(Aest=A_est)

# save data
source("Java-saveData.R")
saveData(A_est, S_est)

detach()