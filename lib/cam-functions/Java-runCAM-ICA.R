# This software is under license BSD 3.
# all parameters (X_mask, K) will be preloaded by Java app.
source("CAM-ICA2dim.R")
source("CAM-ICA3dim.R")

if (K != ncol(X_mask)) {
  cat("Error! Tissues Number doesn't match the data\n")
} else {
  if (K == 2) {                         # 2 dimensions   
    results <- CAM.ICA2dim(X_mask,K)
  } else {                              # 3 or more dimensions
    results <- CAM.ICA3dim(X_mask,K)
  }
}

A <- results[[1]]
S <- results[[2]]
final.result <- list(Aest=A)

# save data
source("Java-saveData.R")
saveData(A, t(S))

detach()
