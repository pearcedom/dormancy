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


CAM.ICA2dim <- function(X, K) {

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
L <- 2

cat("preprocessing ...\n")
#calculate the mean of each class if there are more than one sample

ratio <- X[1,]/X[2,]
geneOrderDecrease <- sort(ratio,decreasing=TRUE,index.return=TRUE)
geneOrderIncrease <- sort(ratio,decreasing=FALSE,index.return=TRUE)

#Initialize the number of OVEDEGs
knum <- 64
loop <- 16

sizeRowOld2 <- 0
knumOld2 <- 0
mixSigOld2 <- matrix(numeric(0), 0,0)

cat('perform ICA ...\n')
for (i in 1:loop){
  #define the step size for every iteration
  if (i<4){
    stepSize <- 64
  } else{
    if (i<8){
      stepSize <- 64  
    } else{
      if (i<11){
        stepSize <- 32
      } else{
        stepSize <- 32
      }
    }
  }

  if (i>1){
      knum <- knum+stepSize
  }
  knumHalf <- knum%/%2
  knumMax <- knumHalf
  knumMin <- knum-knumHalf
  indexMin <- geneOrderDecrease$ix[1:knumMax]
  indexMax <- geneOrderIncrease$ix[1:knumMin]
  index <- c(indexMax,indexMin)
  
  if (i!=1){
    sigMax <- max(X)
  }
  
  nICAResult <- nICA(X[,index])
  icaSig <- unlist(nICAResult[[1]])
  A <- unlist(nICAResult[[2]])
  W <- unlist(nICAResult[[3]])
  sizeRow <- unlist(nICAResult[[4]])
  icaSigAll <- W%*%X
  
  if (i==1){
    AOld <- A
    WOld <- W
    sizeRowOld <- sizeRow
    knumOld <- knum
    sizeRowArray <- sizeRow
    knumArray <- knum
    WArray <- list(W)
    mixSigOld <- X
    indAll <- list(index)
  } else {
    sizeRowOld2 <- sizeRowOld
    AOld2 <- AOld
    WOld2 <- WOld
    knumOld2 <- knumOld
    mixSigOld2 <- mixSigOld
    
    AOld <- A
    WOld <- W
    sizeRowOld <- sizeRow
    knumOld <- knum
    mixSigOld <- X
    sizeRowArray <- c(sizeRowArray,list(sizeRow))
    knumArray <- c(knumArray,knum)
    WArray <- c(WArray,list(W))
    indAll <- c(indAll,list(index))
  }
  if ((sizeRow>30)&(sizeRowOld2!=0)){
    break
  }
}
Wfinal <- WArray[[i-1]]
icaSigfinal <- Wfinal%*%X
icaSigfinal <- Preprocess(icaSig)
A <- solve(Wfinal)

scale <- ginv(A)%*%matrix(1,L,1)
scale <- as.vector(scale)
A <- A%*%diag(scale)

S <- matrix(0,nrow=dim(A)[2],ncol=dim(X)[2])
for (j in 1:ncol(X)){
  S[,j] <- coef(nnls(A,X[,j]))
}

#A <- A/rowSums(A)
#S <- Wfinal %*% X_mask
cat("done.\n")

indexFinal <- indAll[[i-1]]
# windows(width=10,height=10)
plot(icaSigAll[1,],icaSigAll[2,],pch=19,col='blue',cex=0.6,
     main='Scatter Plot of Final Demixing Signal',xlab='S1',ylab='S2',
     xlim=c(0,15),ylim=c(0,15)) 
points(icaSigAll[1,indexFinal],icaSigAll[2,indexFinal],pch=21,col='red',cex=1)

return(list(A,S))
}
