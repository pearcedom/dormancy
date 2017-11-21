# ---------------------------------------------------------------------------------------------
# CM (Compartmental Modeling) algorithm in R version.
# Estimation of kinetic parameters for a two-tissue compartment model.
#
# INPUTS:  X         - M X N matrix where M is the number of image time series, and N is the 
#                      number of pixels; ROI-outlined dynamic imaging data, each column is the 
#                      measured TC curve of a pixel in the ROI. 
#          TCs       - The estimated Tracer/contrast concentrations (TCs) from CAM. each column 
#                      is a compartment TC, for example, when J=3, TCs = [TC_f, TC_s, TC_p].
#          initk     - Stored as an (L-1)x2 matrix, with each row is the initialized 
#                      [Ktrans, Kep], for a particular compartment.
#          del_t     - Time interval between consecutive dynamic images
#                      
# OUTPUTS: eKtrans   - Estimated compartmental volume transfer constants, stored as a 1x(J-1) 
#                      vector with each entry corresponding to a compartment.
#          eKep      - Estimated flux rate constants, stored as a 1x(J-1) vector with each entry 
#                      corresponding to a compartment.
#          pixelwise_Ktran - pixel wise Ktrans.
#
# Author:                         Fan Meng (mengfan@vt.edu)
# Version:                        4.10.2012   
# Original Matlab Version Author: Li Chen (chen2000@vt.edu)
#                                 Tsung-Han Chan (thchan@ieee.org)
# This software is under license BSD 3.
# ---------------------------------------------------------------------------------------------

CM <- function(X, TCs, initk, del_t) { 
  # source the helper functions
  source("functions/nnls_wrapper.R")
  
  L <- nrow(TCs); J <- ncol(TCs)
  Cp <- TCs[,J,drop=FALSE]
  
  # Cost function
  costfun <- function(par) {
    ki <- par[1]; ko <- par[2]
    H  <- matrix(0,L,L)
    for (t in 1:L) {
      H[t,1]   <- ki * exp(-(t-1) * del_t * ko)
      if (t < L) {
        for (tt in 1:(L-t))
          H[t+tt,tt+1] <- H[t,1]
      }
    }
    R <- max(svd(TC - H %*% Cp)$d)
    return(R)
  }
  
  # Optimization
  cat("CM Step 1/1: Optimization\n")
  eKtrans <- rep(0,J-1)
  eKep    <- rep(0,J-1)
  for (i in 1:(J-1)) {
    TC <- TCs[,i]
    x0 <- initk[i,]
    x_optimal  <- optim(x0, costfun, lower=rep(0,2), method="L-BFGS-B")$par
    eKtrans[i] <- x_optimal[1]
    eKep[i]    <- x_optimal[2]
  }
   
  pixelwise_Ktran <- nnls_wrapper(X, cbind(t(t(TCs[,1:(J-1)])/eKtrans), Cp))
  return(list(eKtrans, eKep, pixelwise_Ktran))
}