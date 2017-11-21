# ---------------------------------------------------------------------------------------------
# CAM (Convex Analysis of Mixture) algorithm in R version.
# To find the N corner points of a convex hull formed by the given data X.
#
# INPUTS:  X         - M X N matrix where M is the number of image time series, and N is the 
#                      number of pixels; ROI-outlined dynamic imaging data, each column is the 
#                      measured TC curve of a pixel in the ROI 
#          K         - the number of organs/tissues (or compartments) to be extracted 
#                      (maximally 10)
#          denoise   - whether using multivariate clustering to denoise (denoise=1) or not,
#                      usually for DCEMRI, set denoise=1; for opitcal imaging, set denoise=0;
#          vis       - if using  denoise, whether do visualization of the convexity
#                      
# OUTPUTS: A_est     - M X K matrix; the tracer/contrast concentration of the major organs 
#          S_est     - K X N matrix; estimated spatial distribution maps of the 
#                      organs/tissues/compartments
#
# Author:                         Fan Meng (mengfan@vt.edu) 
#                                 Niya Wang (wangny@vt.edu)
# Version:                        6.4.2013   
# Original Matlab Version Author: Li Chen (chen2000@vt.edu)
#                                 Tsung-Han Chan (thchan@ieee.org)
# This software is under license BSD 3.
# ---------------------------------------------------------------------------------------------

CAM <- function(X,K,denoise,vis) {
  
  # load required "nnls" packages
  # packageExist <- require("nnls")
  # if (!packageExist) {
  #   install.packages("nnls")
  #   library("nnls")
  # }
  
  # load the library "MASS" for the function "ginv"
  packageExist <- require("MASS")
  if (!packageExist) {
    install.packages("MASS")
    library("MASS")
  }
  
  # source the helper functions
  source("functions/affinity.R")
  source("functions/PCA.R")
  source("functions/SL_EM.R")
  source("functions/measure_conv.R")
  source("functions/nnls_wrapper.R")
  
  # use CAM with denoise
  if (denoise == 1) {
    
    real_max <- (.Machine$double.xmax)^0.5  
    real_min <- (.Machine$double.xmin)^0.5
    
    # normalization
    cat("CAM Step 1/4: Normalization\n")
    X_mask <- X
    M <- nrow(X_mask); N<- ncol(X_mask)
    X <- t(t(X_mask)/colSums(X_mask))   
    tmp <- colSums(X_mask)
    tmp_zero <- which(tmp == 0)
    if (length(tmp_zero) != 0)
      X <- X[,-tmp_zero]
    X <- t(X)
    N <- nrow(X); M<- ncol(X)
    
    N0 <- 1500 # if data size is larger than N0, Kmeans will be performed, otherwise APC will be performed.
    
    if (N<N0) {   
      # affinity algorithm
      cat("CAM Step 2/4: Affinity Propagation Clustering\n")
      Xdistance <- matrix(0,N,N)
      for (i in 1:N) {
        temp <- t(t(X) - X[i,])
        temp <- rowSums(temp^2)
        Xdistance[,i] <- temp^0.5
      }
      S <- -Xdistance
      Skk <- min(S)
      for (k in 1:N)
        S[k,k] <- Skk
      
      E <- affinity(X,S,100)
      dE <- diag(E)
      Iexamplar <- which(dE > 0)  # the index of cluster center
      
      R <- 0; eigV <- 0
      # do visualization of the convexity
      if (vis == 1) {
        PCA.result <- PCA(X,2)
        R    <- PCA.result[[1]]
        eigV <- PCA.result[[2]]  
        
        # uncomment the following commad if you want to run this R function alone (without 
        # the use of Java GUI interface), and draw the figure of visualization of the 
        # convexity.
        # windows()
        plot(R[,1], R[,2], xlab="", ylab="", asp=1, col="blue", pch=20)
        examplar <- rep(0,N)
        for (i in 1:N) {
          tmp <- sort(E[i,], decreasing=TRUE, index.return=TRUE)
          examplar[i] <- tmp$ix[1]
        }
        for (i in 1:N)
          lines(c(R[i,1],R[examplar[i],1]), c(R[i,2],R[examplar[i],2]), col="blue")
        lines(R[Iexamplar,1], R[Iexamplar,2])
        lines(R[Iexamplar,1], R[Iexamplar,2], type="p", col="red", pch=1)
      }
      
      
      X <- t(X)
      K0 <- length(Iexamplar)
      init.W <- matrix(1,1,K0)/K0
      init.M <- X[,Iexamplar]
      CV <- cov(t(X))
      init.V <- array(dim=c(dim(CV),K0))
      for (i in 1:K0)
        init.V[,,i] <- CV
      
    } else {
      # kmeans algorithm
      cat("CAM Step 2/4: Kmeans Clustering\n")
      cluster_num <- min(floor(N/50),50)
      
      Kmean_cluster <- kmeans(X,cluster_num,iter.max=100)
      for (i in 1:50){
        tmp <- kmeans(X,cluster_num,iter.max=100)
        if (Kmean_cluster$tot.withinss>tmp$tot.withinss){
          Kmean_cluster <- tmp
        }
      }
      
      
      # do visualization of the convexity
      if (vis == 1) {
        PCA.result <- PCA(X,2)
        R    <- PCA.result[[1]]
        eigV <- PCA.result[[2]]  
        
        # uncomment the following commad if you want to run this R function alone (without 
        # the use of Java GUI interface), and draw the figure of visualization of the 
        # convexity.
        # windows()
        plot(R[,1], R[,2], xlab="", ylab="", asp=1, col="blue", pch=20)
        
        center_index <- c()
        for (i in 1:cluster_num){
          distance <- rowSums((Kmean_cluster$centers[i]-X[which(Kmean_cluster$cluster==i),])^2)
          center_index <- rbind(center_index,which(Kmean_cluster$cluster==i)[which.min(distance)]) 
        }     
        
        for (i in 1:N)
          lines(c(R[i,1],R[center_index[Kmean_cluster$cluster[i]],1]), 
                c(R[i,2],R[center_index[Kmean_cluster$cluster[i]],2]), col="blue")
        
        lines(R[center_index,1], R[center_index,2])
        lines(R[center_index,1], R[center_index,2], type="p", col="red", pch=1)
      }
      
      X <- t(X)
      K0 <- dim(Kmean_cluster$centers)[1]
      init.W <- matrix(1,1,K0)/K0
      init.M <- t(Kmean_cluster$centers)
      CV <- cov(t(X))
      init.V <- array(dim=c(dim(CV),K0))
      for (i in 1:K0)
        init.V[,,i] <- CV
      
      
    }

    
    # EM algorithm
    cat("CAM Step 3/4: Expectation-Maximization (EM)\n")    
    SL_EM.result <- SL_EM(X, init.M, init.V, init.W, 0, 0)
    normindic        <- SL_EM.result[[1]]
    indic            <- SL_EM.result[[2]]
    clustered_center <- SL_EM.result[[3]]
    V                <- SL_EM.result[[4]]
    W                <- SL_EM.result[[5]]
    id_record        <- SL_EM.result[[6]]
    return_flag      <- SL_EM.result[[7]]
    
    filter_TAC <- clustered_center   
    
    # minMargin convex optimization algorithm
    cat("CAM Step 4/4: Convex Analysis of Mixtures (CAM)\n")
    measure_conv.result <- measure_conv(filter_TAC, K) # choose 3 corner points of the convex 
                                                       # hull constructed by filter_TAC
    eA        <- measure_conv.result[[1]]
    cornerind <- measure_conv.result[[2]]
    
    # estimated Cf, Cs, Cp
    eC <- matrix(nrow=M,ncol=K)
    for (k in 1:K){
      outlierremoved <- which(normindic[cornerind[k],] > 0.8)
      tmp_result <- normindic[cornerind[k],outlierremoved,drop=FALSE]
      eC[,k] <- t(tmp_result %*% t(X_mask[,outlierremoved]))/sum(tmp_result)
    }
    A_est <- eC
    S_est <- nnls_wrapper(X_mask, eC)    
    return(list(A_est,S_est))     
  } else {  
    # Use CAM without denoise.    
    # Not used in the project.
    M <- nrow(X); L <- ncol(X)
    
    # principal component analysis for dimension reduction
    cat("CAM Step 1/2: Principal Component Analysis (PCA)\n")
    d <- rowMeans(X)
    U <- X - d
    default <- 19
    tmp <- eigen(U %*% t(U))
    C   <- tmp$vectors; C <- C[,1:default]
    w   <- tmp$values;  D <- diag(w[1:default])
    Xd_t <- t(C) %*% U
    
    # algorithm
    cat("CAM Step 2/2: CAM without Multivariate Clustering (MC)\n")
    A  <- matrix(, default+1, 0)
    Xd <- rbind(Xd_t, matrix(1,1,L))
    P  <- diag(default+1)
    index <- matrix(, 1, 0)
    for (i in 1:(default+1)) {
      ind <- which.max(colSums(abs(P %*% Xd)^2)^0.5)
      A   <- cbind(A, Xd[,ind])
      P   <- diag(default+1) - A %*% ginv(A)
      index <- cbind(index, ind)
    }
    
    A_est <- C %*% Xd_t[,index] + d
    S_est <- nnls_wrapper(Xd_t, Xd_t[,index])
    if (M == 256) {
      index <- c(4:7, 11, (1:(default+1))[-c(4:7,11)])
    } else {
      if (M == 316) {
        index <- c(1:7, 11, 13, 16, (1:(default+1))[-c(1:7,11,13,16)])
      }
    }
    index <- index[1:K]
    
    A_est <- A_est[,index]
    S_est <- S_est[index,]
    return(list(A_est,S_est))
  }
}