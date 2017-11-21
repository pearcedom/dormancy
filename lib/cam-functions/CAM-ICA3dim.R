
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
# Version:                        11/06/2012   
# This software is under license BSD 3.

# ---------------------------------------------------------------------------------------------


CAM.ICA3dim <- function(X, K) {

source('functions2/nICA.R')
source('functions2/Gradient.R')
source('functions2/Whiten.R')
source('functions2/Preprocess.R')
source('functions2/remmean.R')


packageExist <- require("MASS")
if (!packageExist) {
  install.packages("MASS")
  library("MASS")
}

packageExist <- require("nnls")
if (!packageExist) {
  install.packages("nnls")
  library("nnls")
}

X <- t(X)

sizeData <- dim(X)[2]
L <- dim(X)[1]
denom <- as.matrix(colSums(X))
num <- dim(X)[1]
denom <- t(denom[,rep(1,num)])

X_proj <- X/denom

center <- rowMeans(X_proj)

loop <- 48

dis <- sqrt(colSums((X_proj-center)^2))
radius <- max(dis)

sizeIndex <- matrix(numeric(0),0,0)
sizeRowOld2 <- 0
radiusOld2 <- 0


cat('perform ICA ...\n')
for (i in 1:loop){ 
  if (i<2){
    stepSize <- max(dis)/5
  } else{
    if (i<8){
      stepSize <- max(dis)/15 
    } else{
      if (i<12){
        stepSize <- max(dis)/30
      } else{
        if (i<16){
          stepSize <- max(dis)/60
        } else{
          stepSize <- max(dis)/80
        }            
      }
    }
  }
  
  # to reduce the radius of the circle
  radius <- radius-stepSize
  
  index <- matrix(numeric(0),0,0)
  
  dislogical <- dis>=radius
  for (j in 1:length(dis)){
    if (dislogical[j]==TRUE){
      index <- c(index,j)
    }
  }
  

  nICAResult <- nICA(X[,index])
  icaSig <- unlist(nICAResult[[1]])
  A <- unlist(nICAResult[[2]])
  W <- unlist(nICAResult[[3]])
  sizeRow <- unlist(nICAResult[[4]])
  icaSigAll <- W%*%X
  

  
  if (i==1){
    AOld2 <- A
    WOld2 <- W
    sizeRowOld <- sizeRow
    sizeRowArray <- sizeRow
    radiusOld <- radius
    radiusArray <- radius
    WArray <- list(W)
    indAll <- list(index)
    sizeIndex <- c(sizeIndex,length(index))
  } else {
    sizeRowOld2 <- sizeRowOld
    sizeRowOld <- sizeRow
    sizeRowArray <- c(sizeRowArray,sizeRow)
    sizeIndex <- c(sizeIndex,length(index))
    WArray <- c(WArray,list(W))
    AOld <- AOld2
    WOld <- WOld2
    radiusOld <- radiusOld2              
    AOld2 <- A
    WOld2 <- W
    radiusOld2 <- radius
    mixSigOld2 <- X
    radiusArray <- c(radiusArray,radius)
    indAll <- c(indAll,list(index))
  }
    
  if(((abs(sizeRow-sizeRowOld2)>20) | (sizeRow>50))){
    break
  }
}  

Wfinal <- WArray[[i-1]]
A <- solve(W)

scale <- ginv(A)%*%matrix(1,L,1)
scale <- as.vector(scale)
A <- A%*%diag(scale)

S <- matrix(0,nrow=dim(A)[2],ncol=dim(X)[2])
for (j in 1:ncol(X)){
  S[,j] <- coef(nnls(A,X[,j]))
}


return(list(A,S))
}













