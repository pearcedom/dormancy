
# This software is under license BSD 3.
# ------------------ Perform CAM-nWCA Algorithm ---------------------------------------
source("CAM-nWCA.R")
CAM.result <- CAM.nWCA(t(X_mask),K) # all parameters will be preloaded by Java app.
A_est <- CAM.result[[1]]
S_est <- CAM.result[[2]]


# pass results to Java environment
final.result <- list(Aest=A_est)

# save data
source("Java-saveData.R")
saveData(A_est, t(S_est))


detach()
