MDL <- function(X,Xc,S,Sc,A,K){
  
  cat("Calculating MDL")
  L<-dim(X)[1]
  data_size <- dim(X)[2]
  J <- dim(Xc)[2]
  likelihood1 <- (L*J)/2*log(var(as.vector(Xc-A%*%Sc)))
  likelihood2 <- (L*data_size)/2*log(var(as.vector(X-A%*%S)))
  
  sigma1 <- var(as.vector(Xc-A%*%Sc))
  sigma2 <- var(as.vector(X-A%*%S)) 
  
  penalty1 <- (K_est*L)/2*log(J)+(K_est*J)/2*log(L) #+log(J*L)/2
  penalty2 <- (K_est*L)/2*log(data_size)+(K_est*data_size)/2*log(L) #+log(J*L)/2
  
  MDL1 <- likelihood1+penalty1
  MDL2 <- likelihood2+penalty2
  
  return(list(likelihood1,likelihood2,penalty1,penalty2,MDL1,MDL2,sigma1,sigma2))
  
}