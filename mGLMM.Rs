###verification of algorithm###
###K=2, one-dimsional random effects in one response### 
library(statmod)
source("CalculateB.Rs")
##20230625##
###############################################################################

tr <- function(mA) {
  return(sum(diag(mA)))
}

###############################################################################

GVA.mGLMM.f <- function(vtheta,family,vy,mX,id)  
{
  m <- length(unique(id))
  n <- length(which(id==1))
  
  K <- ncol(vy)
  Kpr <- K
  
  p<-0
  for(k in 1:K){
    p <- p + ncol(mX[[k]])
  }
  
  Dinds <- Kpr*((1:Kpr)-1) + (1:Kpr) 	
  Rinds <- which(upper.tri(diag(0,Kpr),diag=TRUE))	 
  Linds <- which(lower.tri(diag(0,Kpr),diag=TRUE))	  
  q <- length(Rinds)
  
  ###########################################################################
  
  vbeta <- matrix(vtheta[1:p],nrow=K)
  vr    <- vtheta[p+(1:q)]
  mR <- matrix(0,Kpr,Kpr)
  mR[Rinds] <- vr
  mR[Dinds] <- exp(mR[Dinds])
  mSigma.inv <- t(mR)%*%mR
  
  sl <- 0
  vmu.X <- matrix(0,nrow=m*n,ncol=K)
  for(k in 1:K){
    vmu.X[,k]   <- mX[[k]]%*%vbeta[,k]
  }
  vmu.Z   <- matrix(0,nrow=m*n,ncol=K)
  vsigma2 <- matrix(0,nrow=m*n,ncol=K)
  mM  <- matrix(0,Kpr,Kpr)
  mRi <- matrix(0,Kpr,Kpr)
  for (i in 1:m) {
    vmui <- vtheta[p + q + (i-1)*(Kpr + Kpr) + (1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    mM <- mM + vmui%*%t(vmui) + t(mRi)%*%mRi
    sl <- sl + sum(log(mRi[Dinds]))
    Lambdai <- t(mRi)%*%mRi
    
    for(k in 1:K){
      vmuiK <- vmui[k] ###matrix(vmui,nrow=2) 
      LambdaiK <- Lambdai[k,k]
      for(j in 1:n){
        index <- (i-1)*n +j
        vmu.Z[index,k] <- vmuiK
        vsigma2[index,k] <- LambdaiK
      }
    }
  }
  veta <- vmu.X + vmu.Z
  ###########################################################################
  family <- "POISSON"
  vB0 <- B0.fun(family,veta,vsigma2) 
  sl <- sl + sum(vy*veta - vB0)
  sl <- sl + m*sum(log(mR[Dinds])) + 0.5*m*Kpr - 0.5*tr(mM%*%mSigma.inv)
  if (family=="POISSON") {
    sl <- sl - sum(lgamma(vy+1))
  }
  return(sl)
}

#############################################################################
GVA.mGLMM.vg <- function(vtheta,family,vy,mX,id) 
{ 
  m <- length(unique(id))
  n <- length(which(id==1))
  
  K <- ncol(vy)
  Kpr <- K
  
  p<-0
  for(k in 1:K){
    p <- p + ncol(mX[[k]])
  }
  
  Dinds <- Kpr*((1:Kpr)-1) + (1:Kpr) 	
  Rinds <- which(upper.tri(diag(0,Kpr),diag=TRUE))	 
  Linds <- which(lower.tri(diag(0,Kpr),diag=TRUE))	  
  q <- length(Rinds)
  
  ###########################################################################
  
  vbeta <- matrix(vtheta[1:p],nrow=K)
  vr    <- vtheta[p+(1:q)]
  mR <- matrix(0,Kpr,Kpr)
  mR[Rinds] <- vr
  mR[Dinds] <- exp(mR[Dinds])
  mSigma.inv <- t(mR)%*%mR
  
  ###########################################################################
  vmu.X <- matrix(0,nrow=m*n,ncol=K)
  for(k in 1:K){
    vmu.X[,k]   <- mX[[k]]%*%vbeta[,k]
  }
  
  vmu.Z   <- matrix(0,nrow=m*n,ncol=K)
  vsigma2 <- matrix(0,nrow=m*n,ncol=K)
  
  mM  <- matrix(0,Kpr,Kpr)
  mRi <- matrix(0,Kpr,Kpr)
  for (i in 1:m) {
    vmui <- vtheta[p + q + (i-1)*(Kpr + Kpr) + (1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    mM <- mM + vmui%*%t(vmui) + t(mRi)%*%mRi
    Lambdai <- t(mRi)%*%mRi
    
    for(k in 1:K){
      vmuiK    <- vmui[k] 
      LambdaiK <- Lambdai[k,k]
      
      for(j in 1:n){
        index <- (i-1)*n +j
        vmu.Z[index,k] <- vmuiK
        vsigma2[index,k] <- LambdaiK
      }
    }
  }
  veta <- vmu.X + vmu.Z
  
  ###########################################################################    
  
  lvB  <- B1234.fun(family,veta,vsigma2)
  vB1 <- lvB[[1]]
  vB2 <- lvB[[2]] #; dim(vB2) <- NULL
  vB3 <- lvB[[3]] #; dim(vB3) <- NULL
  vB4 <- lvB[[4]] #; dim(vB4) <- NULL
  veps <- vy - vB1 
  
  ###########################################################################
  vg.vbeta <- c() 
  
  for(k in 1:K){
    vg.vbeta <- c(vg.vbeta,apply(mX[[k]]*veps[,k],2,sum))
  }
  
  r11 <- log(mR[1,1])
  r12 <- mR[1,2]
  r22 <- log(mR[2,2])
  vg.r11 <- 0
  vg.r12 <- 0
  vg.r22 <- 0
  
  for(i in 1:m){
    vmui <- vtheta[p+q+(i-1)*(Kpr+Kpr)+(1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]
    
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    
    vmui1 <- vmui[1]
    vmui2 <- vmui[2]
    ri11  <- log(mRi[1,1])
    ri22  <- log(mRi[2,2])
    
    vg.r11 <- vg.r11 + 1 - exp(2*r11)*vmui1^2 - r12*exp(r11)*vmui1*vmui2 - exp(2*r11+2*ri11)
    vg.r12 <- vg.r12 - exp(r11)*vmui1*vmui2   - r12*vmui2^2 - r12*exp(2*ri22)
    vg.r22 <- vg.r22 + 1 - exp(2*r22)*vmui2^2 - exp(2*r22)*exp(2*ri22)
  }
  
  vg.mR <- c(vg.r11,vg.r12,vg.r22)
  vg.vtheta <- c(vg.vbeta,vg.mR)
  
  ###########################################################################
  vg <- vg.vtheta
  count <- 0
  linds <- list()
  for (i in 1:m) {
    linds[[i]] <- count + (1:n)
    count <- count + n
  }
  
  lvg.vxi <- list()
  for(i in 1:m){
    inds <- linds[[i]]
    n   <- length(inds)
    
    vmui <- vtheta[p+q+(i-1)*(Kpr+Kpr)+(1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]
    
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    
    ri11 <- log(mRi[1,1])
    ri22 <- log(mRi[2,2])
    
    mXi1 <- mX[[1]][inds,]
    vB1i1 <- vB1[inds,1]
    vB2i1 <- vB2[inds,1]
    vB3i1 <- vB3[inds,1]
    vB4i1 <- vB4[inds,1]
    
    mXi2 <- mX[[2]][inds,]
    vB1i2 <- vB1[inds,2]
    vB2i2 <- vB2[inds,2]
    vB3i2 <- vB3[inds,2]
    vB4i2 <- vB4[inds,2]
    
    vg.vmui1 <- sum(vy[inds,1] - vB1i1) - vmui1*exp(2*r11)   - r12*vmui2*exp(r11) 
    vg.vmui2 <- sum(vy[inds,2] - vB1i2) - r12*vmui1*exp(r11) - (r12^2+exp(2*r22))*vmui2
    
    vg.ri11 <- -sum(vB2i1)*exp(2*ri11) + 1 - exp(2*r11 + 2*ri11)
    vg.ri22 <- -sum(vB2i2)*exp(2*ri22) + 1 - (r12^2 + exp(2*r22))*exp(2*ri22)
    
    vg.xi <- c(vg.vmui1,vg.vmui2,vg.ri11,vg.ri22)
    
    lvg.vxi[[i]] <- matrix(vg.xi)
    vg <- c(vg,lvg.vxi[[i]])
  }
  
  res.vg <- list(vg.0=vg.vtheta,lvg.i=lvg.vxi,vg=vg)
  
  return(res.vg)
}

###Please remember to start here###
###############################################################################
GVA.grouped.vgANDmH <- function(vtheta,family,vy,mX,id) 
{ 
  m <- length(unique(id))
  n <- length(which(id==1))
  
  K <- length(mX)
  Kpr <- K
  
  p<-0
  for(k in 1:K){
    p <- p + ncol(mX[[k]])
  }
  
  Dinds <- Kpr*((1:Kpr)-1) + (1:Kpr) 	
  Rinds <- which(upper.tri(diag(0,Kpr),diag=TRUE))	 
  Linds <- which(lower.tri(diag(0,Kpr),diag=TRUE))	  
  q <- length(Rinds)
  
  ###########################################################################
  vbeta <- matrix(vtheta[1:p],nrow=2)
  vr    <- vtheta[p+(1:q)]
  mR    <- matrix(0,Kpr,Kpr)
  mR[Rinds]  <- vr
  mR[Dinds]  <- exp(mR[Dinds])
  mSigma.inv <- t(mR)%*%mR
  
  ###########################################################################
  
  vmu.X <- matrix(0,nrow=m*n,ncol=K)
  for(k in 1:K){
    vmu.X[,k]   <- mX[[k]]%*%vbeta[,k]
  }
  
  vmu.Z   <- matrix(0,nrow=m*n,ncol=K)
  vsigma2 <- matrix(0,nrow=m*n,ncol=K)
  mM  <- matrix(0,Kpr,Kpr)
  mRi <- matrix(0,Kpr,Kpr)
  
  
  for (i in 1:m) {
    vmui <- vtheta[p + q + (i-1)*(Kpr + Kpr) + (1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]
    
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    mM <- mM + vmui%*%t(vmui) + t(mRi)%*%mRi
    
    Lambdai <- t(mRi)%*%mRi
    for(k in 1:K){
      vmuiK <- vmui[k]
      LambdaiK <- Lambdai[k,k]
      
      for(j in 1:n){
        index <- (i-1)*n +j
        vmu.Z[index,k] <- vmuiK
        vsigma2[index,k] <- LambdaiK
      }
    }
  }
  veta <- vmu.X + vmu.Z
  
  ###########################################################################    
  lvB  <- B1234.fun(family,veta,vsigma2)
  vB1 <- lvB[[1]]
  vB2 <- lvB[[2]] #; dim(vB2) <- NULL
  vB3 <- lvB[[3]] #; dim(vB3) <- NULL
  vB4 <- lvB[[4]] #; dim(vB4) <- NULL
  veps <- vy - vB1 
  
  ###########################################################################
  res.vg <- GVA.mGLMM.vg(vtheta,family,vy,mX,id) 
  vg.vtheta <- res.vg$vg.0
  lvg.vxi   <- res.vg$lvg.i 
  ###########################################################################
  mH.beta.beta <- matrix(0,4,4)
  
  mH.beta.beta[1:2,1:2] <-  -t(mX[[1]]*vB2[,1])%*%mX[[1]]
  mH.beta.beta[3:4,3:4] <-  -t(mX[[2]]*vB2[,2])%*%mX[[2]]
  
  mH.mR.mR <- matrix(0,q,q)
  r11 <- log(mR[1,1])
  r12 <- mR[1,2]
  r22 <- log(mR[2,2])
  mH.r11.r11 <- 0
  mH.r11.r12 <- 0
  mH.r11.r22 <- 0
  mH.r12.r12 <- 0
  mH.r12.r22 <- 0
  mH.r22.r22 <- 0
  
  for(i in 1:m){
    vmui <- vtheta[p + q + (i-1)*(Kpr + Kpr) +(1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]

    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    
    ri11 <- log(mRi[1,1])
    ri22 <- log(mRi[2,2])
    
    vmui1 <- vmui[1] 
    vmui2 <- vmui[2] 
    
    mH.r11.r11 <- mH.r11.r11 - 2*exp(2*r11)*vmui1^2 - r12*exp(r11)*vmui1*vmui2 - 2*exp(2*r11+2*ri11)
    mH.r11.r12 <- mH.r11.r12 - exp(r11)*vmui1*vmui2
    mH.r12.r12 <- mH.r12.r12 - vmui2^2 - exp(2*ri22) 
    mH.r22.r22 <- mH.r22.r22 - 2*vmui2^2*exp(2*r22) - 2*exp(2*r22)*exp(2*ri22)
  }
  
  mH.mR.mR[1,1] <- mH.r11.r11
  mH.mR.mR[1,2] <- mH.r11.r12
  mH.mR.mR[1,3] <- mH.r11.r22
  mH.mR.mR[2,2] <- mH.r12.r12
  mH.mR.mR[2,3] <- mH.r12.r22
  mH.mR.mR[3,3] <- mH.r22.r22
  
  mH.mR.mR[2,1] <- mH.mR.mR[1,2]
  mH.mR.mR[3,1] <- mH.mR.mR[1,3]
  mH.mR.mR[3,2] <- mH.mR.mR[2,3]
  
  ##beta.mR###
  mH.theta.theta <- matrix(0,4+q,4+q)
  mH.theta.theta[1:4,1:4]<- mH.beta.beta
  mH.theta.theta[4+1:q,4+1:q] <- mH.mR.mR
  
  ###mR.xi.xi###
  count <- 0
  linds <- list()
  for (i in 1:m) {
    linds[[i]] <- count + (1:n)
    count <- count + n
  }
  
  
  lmH.ii <- list()
  lmH.0i <- list()
  
  for(i in 1:m){
    
    inds <- linds[[i]]
    ni   <- length(inds)
    
    vmui <- vtheta[p+q+(i-1)*(Kpr+Kpr)+(1:Kpr)]
    vri  <- vtheta[p + q + (i-1)*(Kpr + Kpr) + Kpr + (1:Kpr)]
     
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    
    ri11 <- log(mRi[1,1])
    ri22 <- log(mRi[2,2])
    
    mXi1 <- mX[[1]][inds,]
    vB1i1 <- vB1[inds,1]
    vB2i1 <- vB2[inds,1]
    vB3i1 <- vB3[inds,1]
    vB4i1 <- vB4[inds,1]
    
    mXi2 <- mX[[2]][inds,]
    vB1i2 <- vB1[inds,2]
    vB2i2 <- vB2[inds,2]
    vB3i2 <- vB3[inds,2]
    vB4i2 <- vB4[inds,2]
    
    ###mH.mRi.mRi###
    mH.mRi.mRi <- matrix(0,Kpr,Kpr)
    
    mH.ri11.ri11 <- -sum(vB4i1)*exp(4*ri11) - sum(vB2i1)*2*exp(2*ri11) - 2*exp(2*r11+2*ri11)
    mH.ri11.ri22 <- 0
    mH.ri22.ri22 <- -sum(vB4i2)*exp(4*ri22) - sum(vB2i2)*2*exp(2*ri22) - 2*(r12^2+exp(2*r22))*exp(2*ri22)
    
    mH.mRi.mRi[1,1] <- mH.ri11.ri11
    mH.mRi.mRi[2,2] <- mH.ri22.ri22
    
    ###mH.vmui.vmui###
    mH.vmui.vmui <- matrix(0,2,2)
    mH.vmui.vmui[1,1] <- -sum(vB2i1) - exp(2*r11)
    mH.vmui.vmui[1,2] <- -r12*exp(r11)
    mH.vmui.vmui[2,1] <- mH.vmui.vmui[1,2] 
    mH.vmui.vmui[2,2] <- -sum(vB2i2) - r12^2 - exp(2*r22)
    
    ###mH.vmui.mRi###
    mH.vmui.mRi <- matrix(0,2,Kpr)
    mH.vmui.mRi[1,1] <- -sum(vB3i1)*exp(2*ri11)
    mH.vmui.mRi[2,2] <- -sum(vB3i2)*exp(2*ri22)
    
    mH.xi.xi <- matrix(0, 2 + Kpr, 2 + Kpr)
    mH.xi.xi[1:2,1:2]   <- mH.vmui.vmui
    mH.xi.xi[1:2,2+1:Kpr] <- mH.vmui.mRi
    mH.xi.xi[2+1:Kpr,1:2] <- t(mH.xi.xi[1:2,2+1:Kpr])
    mH.xi.xi[2+1:Kpr,2+1:Kpr] <- mH.mRi.mRi
    
    ##mH.theta.xi###
    mH.theta.xi <- matrix(0, 4 + q, 2 + Kpr)
    mH.beta1.vmui1 <- -apply(mXi1*vB2i1,2,sum)
    mH.beta1.ri11  <- -apply(mXi1*vB3i1,2,sum)*exp(2*ri11)
    
    mH.beta2.vmui2 <- -apply(mXi2*vB2i2,2,sum)  
    mH.beta2.ri22  <- -apply(mXi2*vB3i2,2,sum)*exp(2*ri22)  
    
    mH.beta.xi <- matrix(0, 4, 2+Kpr)
    mH.beta.xi[1:2,1] <- mH.beta1.vmui1
    mH.beta.xi[1:2,3] <- mH.beta1.ri11
    
    mH.beta.xi[3:4,2] <- mH.beta2.vmui2
    mH.beta.xi[3:4,4] <- mH.beta2.ri22
    
    mH.mR.xi <- matrix(0,q,2+Kpr)
    mH.r11.vmui1 <- -2*vmui1*exp(2*r11) - r12*exp(r11)*vmui2
    mH.r11.vmui2 <- -r12*exp(r11)*vmui1
    mH.r11.ri11  <- -2*exp(2*r11+2*ri11)
    mH.r11.ri22  <-  0
    
    mH.r12.vmui1 <- -vmui2*exp(r11)
    mH.r12.vmui2 <- -vmui1*exp(r11) - 2*r12*vmui2
    mH.r12.ri11  <- 0
    mH.r12.ri22  <- -2*r12*exp(2*ri22)
    
    mH.r22.vmui1 <- 0
    mH.r22.vmui2 <- -2*vmui2*exp(2*r22)
    mH.r22.ri11  <- 0
    mH.r22.ri22  <- -2*exp(2*r22)*exp(2*ri22)
    
    mH.mR.xi[1,] <- c(mH.r11.vmui1,mH.r11.vmui2,mH.r11.ri11,mH.r11.ri22)
    mH.mR.xi[2,] <- c(mH.r12.vmui1,mH.r12.vmui2,mH.r12.ri11,mH.r12.ri22)
    mH.mR.xi[3,] <- c(mH.r22.vmui1,mH.r22.vmui2,mH.r22.ri11,mH.r22.ri22)
    
    mH.theta.xi[1:4,1:(2+Kpr)] <- mH.beta.xi
    mH.theta.xi[4+1:q,1:(2+Kpr)] <- mH.mR.xi
    
    lmH.ii[[i]] <- mH.xi.xi
    lmH.0i[[i]] <- mH.theta.xi
  }
  
  res.mH <- list(mH.00=mH.theta.theta,lmH.0i=lmH.0i,lmH.ii=lmH.ii)
  
  return(list(res.vg=res.vg,res.mH=res.mH))
}


###############################################################################

GVA.mGLMM.NRDIR <- function(vtheta,family,vy,mX,id)    
{ 
  
  m <- length(unique(id))
  n <- length(which(id==1))
  
  K <- length(mX)
  Kpr <- K
  
  res <- GVA.grouped.vgANDmH(vtheta,family,vy,mX,id)    
  res.vg <- res$res.vg
  res.mH <- res$res.mH
  vg <- c(res.vg$vg.0)
  vg.til <- res.vg$vg.0
  mH.til <- res.mH$mH.00
  for (i in 1:m) {
    vg.i <- res.vg$lvg.i[[i]]
    mH.0i     <- res.mH$lmH.0i[[i]]
    mH.ii.inv <- solve(res.mH$lmH.ii[[i]],tol=1.0E-99)
    mA <- mH.0i%*%mH.ii.inv
    vg <- c(vg,vg.i)
    vg.til <- vg.til - mA%*%vg.i
    mH.til <- mH.til - mA%*%t(mH.0i)
  } 
  
  sI <- solve(mH.til)
  vs <- sI%*%vg.til
  
  dvtheta <- vs
  for (i in 1:m) {
    dvtheta <- c(dvtheta,solve(res.mH$lmH.ii[[i]],res.vg$lvg.i[[i]] - t(res.mH$lmH.0i[[i]])%*%vs))
  }
  return(list(dvtheta=-dvtheta,sI=-sI,vg=vg))
}

###############################################################################

getTheta <- function(vbeta,mSigma,mu,Lambda) 
{    
  m <- length(mu)
  K <- ncol(mSigma)
  
  Dinds <- K*((1:K)-1) + (1:K) 	
  Rinds <- which(upper.tri(diag(0,K),diag=TRUE))      
  
  mSigma.inv <- solve(mSigma,tol=1.0E-99)  
  mR <- chol(mSigma.inv)
  mR[Dinds] <- log(mR[Dinds])
  
  vtheta <- c(vbeta,mR[Rinds])
  for (i in 1:m) {
    mR <- chol(Lambda[[i]])    
    mR[Dinds] <- log(mR[Dinds])
    vtheta <- c(vtheta,mu[[i]],mR[Dinds])
  }
  return(vtheta)
}

###############################################################################

thetaToPars <- function(vtheta,vy,mX,mu) 
{
  m <- length(mu)
  K <- ncol(mX[[1]])
  
  p<-0
  for(k in 1:K){
    p <- p + ncol(mX[[k]])
  }
  
  Dinds <- K*((1:K)-1) + (1:K) 	
  Rinds <- which(upper.tri(diag(0,K),diag=TRUE))	  
  Linds <- which(lower.tri(diag(0,K),diag=TRUE))	  
  q <- length(Rinds)
  
  vbeta <- vtheta[1:p]
  vr    <- vtheta[p+(1:q)]
  mR <- matrix(0,K,K)
  mR[Rinds] <- vr
  mR[Dinds] <- exp(mR[Dinds])
  mSigma.inv <- t(mR)%*%mR   
  mSigma     <- solve(mSigma.inv,tol=1.0E-99) 
  
  mRi <- matrix(0,K,K)
  lvmu     <- list()
  lmLambda <- list()
  for (i in 1:m) {
    lvmu[[i]] <- vtheta[p+q+(i-1)*(K+K)+(1:K)]
    vri  <- vtheta[p+q+(i-1)*(K+K)+K+(1:K)]
    mRi[Dinds] <- vri
    mRi[Dinds] <- exp(mRi[Dinds])
    lmLambda[[i]] <- t(mRi)%*%mRi
  }
  
  return(list(vbeta=vbeta,mSigma=mSigma,lvmu=lvmu,lmLambda=lmLambda))
}

###############################################################################

GVA.mGLMM.FIT <- function(vbeta,mSigma,vmu,vlambda,family,vy,mX,id) 
{
  MAXITER  <- 100
  EPS.TERM.PAR  <- 1.0E-4
  EPS.TERM.GRAD <- 1.0E-8
  
  family <- "POISSON"
  vy <- simy1$y
  vbeta <- t(vbeta)
  vtheta <- getTheta(vbeta,mSigma,vmu,vlambda) 
  length(vtheta)
  sl <- GVA.mGLMM.f(vtheta,family,vy,mX,id)  
  
  cat("Initial variational log-likelihood = ",sl)
  cat("\n\n")
  cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
  
  for (ITER in 1:MAXITER) {
    #ITER <- 1
    vtheta.old <- vtheta
    
    # Calculate Newton-Raphson Direction 
    res <- GVA.mGLMM.NRDIR(vtheta,family,vy,mX,id)    
    
    stepsize <- 1                
    vtheta.temp <- vtheta + stepsize*res$dvtheta
    sl.temp <- GVA.mGLMM.f(vtheta.temp,family,vy,mX,id)  
    
    if (is.nan(sl.temp)) { 
      DOGRIDSEARCH <- TRUE
    } else {            
      if (sl.temp>sl) { DOGRIDSEARCH <- FALSE } 
      else { DOGRIDSEARCH <- TRUE }
    }         
    
    if (DOGRIDSEARCH) {
      NG <- 15
      vsl <- c()
      
      # Allow for small negative step sizes 
      # (just in case Hessian is negitive definite)
      vstepsize <- seq(-0.25,1,,NG)[-NG]  
      
      for (JTER in 1:length(vstepsize)) {
        vtheta.temp <- vtheta + vstepsize[JTER]*res$dvtheta
        vsl[JTER]   <- GVA.mGLMM.f(vtheta.temp,family,vy,mX,id)  
      }                   
      sl <- max(vsl)
      stepsize <- vstepsize[which.max(vsl)]
      vtheta <- vtheta + stepsize*res$dvtheta
    } else {
      sl <- sl.temp
      vtheta <- vtheta.temp    
    }
    
    # Check termination condtions
    err.par <- max(abs(vtheta[1:7]-vtheta.old[1:7]))/max(abs(vtheta[1:7]))
    err.vg  <- sum(res$vg^2)/length(vtheta)
    
    cat(sprintf("%4d    %E    %+E    %E  %E\n",ITER,sl,stepsize,err.vg,err.par))
    
    if (err.vg<EPS.TERM.GRAD) { break; }
    if (err.par<EPS.TERM.PAR) { break; }
    
    #if (exp(vtheta[5])<5.0E-4) { break; } # Break if sigma2 is getting too small
    #if (exp(vtheta[6])<5.0E-4) { break; } # Break if sigma2 is getting too small
    #if (exp(vtheta[7])<5.0E-4) { break; } # Break if sigma2 is getting too small
  }
  
  pars <- thetaToPars(vtheta,vy,mX,vmu)
  #vbeta.serr <- sqrt(diag(res$sI[1:4,1:4]))
  
  return(list(vbeta=pars$vbeta,mSigma=pars$mSigma,lvmu=pars$lvmu,lmLambda=pars$lmLambda,sl=sl,ITER=ITER))
}

################################################################################




