
# ---------------------------------------------------------------------------------------------
# 
# INPUTS:  X         - M X N matrix where M is the number of image time series, and N is the 
#                      number of pixels; ROI-outlined dynamic imaging data, each column is the 
#                      measured TC curve of a pixel in the ROI 
#          K         - the number of organs/tissues (or compartments) to be extracted 
#                      (maximally 10)
#         
# OUTPUTS: A_est     - M X K matrix; the tracer/contrast concentration of the major organs 
#          S_est     - K X N matrix; estimated spatial distribution maps of the 
#                      organs/tissues/compartments
#
# Author:                         Niya Wang (wangny@vt.edu)
# Version:                        11/02/2012  
# This software is under license BSD 3.
# ---------------------------------------------------------------------------------------------


CAM.nWCA <- function(X, K) {

  # load the library "MASS" for the function "ginv"
  packageExist <- require("MASS")
  if (!packageExist) {
    install.packages("MASS")
    library("MASS")
  }

  # load the library "geometry" for the function "convhulln"
  packageExist <- require("geometry")
  if (!packageExist) {
    install.packages("geometry")
    library("geometry")
  }
  
  # load the library "nnls" for the function "nnls"
  packageExist <- require("nnls")
  if (!packageExist) {
    install.packages("nnls")
    library("nnls")
  }
  
  source("functions3/measure_conv.R")
 
  data_size<- dim(X)[2]
  L <- dim(X)[1]
  
  ##### use kmeans to cluster the observation #####
  cat("\nPerforming kmeans to cluster data ... \n")
  
  denom <- as.matrix(colSums(X))
  num <- dim(X)[1]
  denom <- t(denom[,rep(1,num)])
  
  X_proj <- X/denom
  cluster_num <- 50
  
  cluster <- kmeans(t(X_proj),cluster_num,iter.max=100)
  for (i in 1:50){
    tmp <- kmeans(t(X_proj),cluster_num,iter.max=100)
    if (cluster$tot.withinss>tmp$tot.withinss){
      cluster <- tmp
    }
  }
  
  small_cluster <- matrix(numeric(0),0,0)
  for (k in 1:cluster_num){
    if (cluster$size[k]<0.1*dim(X)[2]/cluster_num)
      small_cluster <- c(small_cluster,k)
  }
  
  if (length(small_cluster)==0){
    cluster <- cluster
  } else {
    cluster$centers <- cluster$centers[-small_cluster,]
  }
  
  J <- dim(cluster$centers)[1]
  
  convex <- convhulln(rbind(cluster$centers,0))
  
  corner <- matrix(numeric(0), 0,0)
  for (i in 1:L){
    corner <- union(corner,convex[,i])
  }
  
  for (j in 1:length(corner)){
    if (corner[j]==(J+1)){
      break
    }
  }
  corner <- corner[-j]  # throw away the origin point
  #corner <- c(1:35)
  J_out <- length(corner)
  
  ##### estimate A and S ###########
  cat("\nEstimating A and S ... \n")  

  
  if (K==J_out){
   A_est <- t(cluster$centers[corner,])
  } else {
    cornerResult <- measure_conv(t(cluster$centers[corner,]),K)
  
    A_est <- cornerResult[[1]]
    ind <- cornerResult[[2]]    
  }
  
  scale <- ginv(A_est)%*%matrix(1,L,1)
  scale <- as.vector(scale)
  A_est <- A_est%*%diag(scale)
  
  S_est <- matrix(0,nrow=dim(A_est)[2],ncol=dim(X_proj)[2])
  for (i in 1:ncol(X_proj)){
  S_est[,i] <- coef(nnls(A_est,X[,i]))
  }
  
  return(list(A_est,S_est))

}



