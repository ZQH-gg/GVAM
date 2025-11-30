rm(list = ls())
getwd()
setwd("C:/Users/admin/Desktop")
library(scales)
library(ggplot2)
library(Matrix)
library(MCMCglmm)
library(mvtnorm) 
source("CalculateB.Rs")
source("GVArandint.Rs")
source("mGLMM.Rs")
################ functions ##############
gendat.mglmm <- function(id, fixed.MM, ran.MM, beta, Sigma, phi = NULL, cutoffs = NULL, family) {
  
  if(!(family %in% c("GAUSSIAN","POISSON","BINARY","ORDINAL"))) 
    stop("Specified family not permitted. Sorry!")
  K <- length(beta)
  for(k in 1:K){
    
    if(!all(fixed.MM[[k]][,1] == 1)) 
      warning("Is your first column of fixed.MM an intercept? Should it be?")
    if(!all(ran.MM[[k]][,1] == 1)) 
      warning("Is your first column of ran.MM an intercept? Should it be?")
    
    if(length(beta[[k]]) != ncol(fixed.MM[[k]])) 
      stop("Dimensions of beta and fixed.MM are inconsistent. Please rectify. Thanks!")
  }
  
  randim <- sum(sapply(ran.MM, ncol))
  if(ncol(Sigma) != randim) 
    stop("Dimensions of Sigma, ran.MM are inconsistent. Please rectify. Thanks!")

  m <- length(unique(id))  
  u <- rmvnorm(m, rep(0,randim), sigma = Sigma+diag(1e-8, ncol(Sigma))) 
  u[,which(diag(Sigma) == 0)] <- 0
 
  group_vector <- rep(1:K, times = sapply(ran.MM, ncol))  
  u_list <- lapply(split(1:ncol(u), group_vector), function(col_indices) {
    u[, col_indices, drop = FALSE] 
  })
  
  sim.y <- matrix(NA, nrow = length(id), ncol = K)
  for(i in 1:m) {
    sel.i <- which(id == i)
    n_i <- length(sel.i)
    for(k in 1:K){
      eta <- fixed.MM[[k]][sel.i, , drop = FALSE] %*% beta[[k]] + ran.MM[[k]][sel.i, , drop = FALSE] %*% u_list[[k]][i,]
      
      if(family == "GAUSSIAN") sim.y[sel.i,k] <- rnorm(n_i, mean = eta, sd = sqrt(phi[k]))
      
      if(family == "POISSON" ) sim.y[sel.i,k] <- rpois(n_i, exp(eta))
      
      if(family == "BINARY"  ) sim.y[sel.i,k] <- rbinom(n_i, 1, 1/(1+exp(-eta))) 
      if(family == "ORDINAL") sim.y[sel.i,k] <- rorddata(cutoffs = cutoffs[k,], eta = eta)
    }
  }
  out <- list(y=sim.y, u=u_list)
  return(out)
}

GVAM.FIT <- function(Xdependent, id, Y, fixed.MM, beta0, gamma0, mu0, zeta0, family) {
  K<-ncol(Y)
  p <- lapply(fixed.MM, ncol)
  
  for(k in 1:K){
    if(!all(fixed.MM[[k]][,1] == 1)) 
      warning("Is your first column of fixed.MM an intercept? Should it be?")
  }
  converged <- rep(FALSE,K)
  times <- rep(0,K+1)
  beta <- lapply(p, function(x) rep(0, x))
  sigma <- rep(0,K)
  combinations <- combn(K, 2, simplify = FALSE)
  rho <- rep(0,length(combinations))

  unique.ids <- unique(id)
  m  <- length(unique.ids)
  gamma <- gamma0 
  vmu   <- matrix(mu0,m,1)
  vzeta <- matrix(zeta0,m,1)
  for (k in 1:K){
    vy <- Y[,k]
    mX <- fixed.MM[[k]]	
    vbeta <- matrix(beta0,p[[k]],1)
    
    cat("Fitting model using GVA", "  k=", k, "\n")						
    bval  <- proc.time()
    res <- try(GVA.fit <- GVA.randint.FIT(id,vy,mX,vbeta,sigma2=exp(gamma),vmu,vlambda=exp(vzeta),family), silent = TRUE) 
    eval  <- proc.time()
    
    if((class(res)[1])!="try-error"){ 
      converged[k] <- TRUE
      times[k] <- eval[3] - bval[3]
      beta[[k]] <- res$vbeta
      sigma[k] <- res$sigma2
    } else {
      converged[k] <- FALSE
    }
  }
  
  
  if (all(converged)) {
    if (Xdependent==TRUE){
      bval  <- proc.time() 
      rho <- sapply(combinations, function(x){
        vcov <- cov(Y[,c(x[1],x[2])])[1,2]
        ve <- exp((sigma[x[1]]+sigma[x[2]])/2)
        vc <- mean(exp(fixed.MM[[x[1]]]%*%beta[[x[1]]] + fixed.MM[[x[2]]]%*%beta[[x[2]]]))
        vd <- mean(exp(fixed.MM[[x[1]]]%*%beta[[x[1]]])) * mean(exp(fixed.MM[[x[2]]]%*%beta[[x[2]]]))
        log((vcov+ve*vd)/(ve*vc))/sqrt(sigma[x[1]]*sigma[x[2]])
      })
      eval  <- proc.time()
    }else {
      bval  <- proc.time()
      rho <- sapply(combinations, function(x){
        vcov <- cov(Y[,c(x[1],x[2])])[1,2]
        ve <- exp((sigma[x[1]]+sigma[x[2]])/2)
        vc <- mean(exp(fixed.MM[[x[1]]]%*%beta[[x[1]]] + fixed.MM[[x[2]]]%*%beta[[x[2]]]))
        log(1+vcov/(ve*vc))/sqrt(sigma[x[1]]*sigma[x[2]])
      })
      eval  <- proc.time()
    }
    names(rho) <- sapply(combinations, function(x) paste0("rho", x[1], x[2]))
    times[K+1] <- eval[3] - bval[3]
  }
  out <- list(converged=converged, time=sum(times), beta=beta, sigma=sigma, rho=rho)
  return(out)
}

################ simulation setting #####################
Xdependent=FALSE 
M <- c(10,50,100,500,1000,3000,5000,9000) 
ID <- lapply(M, function(m) rep(1:m, round(rnorm(m, mean = 10, sd = 1))))
MS <- 100 
K = 2 
P = c(4, 2) 
true.beta <- list(
  beta1 = c(-1.0, 0.5, -0.8, 0.2, 0.6), 
  beta2 = c(-0.5, 0.3, -1.0)  
  )  
true.Sigma <- matrix(c(1,0.5,0.5,0.8),ncol=2)
true.Rho <- true.Sigma/sqrt(diag(true.Sigma) %o% diag(true.Sigma))

################ begin simulating #####################
res.gvam <- array(list(), dim = c(length(M), MS))
res.mcmc <- array(list(), dim = c(length(M), MS))
set.seed(0)
ITER0 <- 1
for (ITER in ITER0:length(M)) {
  m <- M[ITER]
  id <- ID[[ITER]] 
  trial0 <- 1
  for (trial in trial0:MS) {
    
    X <- list()
    for (k in 1:K){
        X[[k]] <- rmvnorm(length(id),rep(0,P[k]),sigma=diag(P[k]))
    }
    
    fixed.MM <- lapply(X, function(mat) {
      cbind(1, mat)
    })
    
    ran.MM <- replicate(K, matrix(rep(1,length(id)),ncol=1), simplify = FALSE)
    simy <- gendat.mglmm(id = id, fixed.MM = fixed.MM, ran.MM = ran.MM, beta = true.beta, Sigma = true.Sigma, family = "POISSON")
    Y <- simy$y  
    
    Xdependent <- Xdependent
    beta0 <- 0  
    gamma0 <- 0  
    mu0   <- 0  
    zeta0 <- 0 
    family = "POISSON"
    cat("m=",m,"    trial=",trial, "\n")
    res.gvam[[ITER,trial]] <- GVAM.FIT(Xdependent=Xdependent, id=id, Y=Y, fixed.MM=fixed.MM, beta0=beta0, gamma0=gamma0, mu0=mu0, zeta0=zeta0, family = family)
    
   
    converged.mcmc <- FALSE
    time.mcmc <- 0
    beta.mcmc <- list(rep(0,5), rep(0,3))
    sigma.mcmc <- rep(0,2)
    rho.mcmc <- 0
    
    df.mcmc<-data.frame(
      vy1 = Y[,1],
      mX11 = fixed.MM[[1]][,2],
      mX12 = fixed.MM[[1]][,3],
      mX13 = fixed.MM[[1]][,4],
      mX14 = fixed.MM[[1]][,5],
      vy2 = Y[,2],
      mX21 = fixed.MM[[2]][,2],
      mX22 = fixed.MM[[2]][,3],
      factor = factor(id),
      ID = c(1:length(id))
    )
    prior = list(
      R = list(V = diag(2), nu = 0.002), 
      G = list(G1 = list(V = diag(2), nu = 0.002))
    )
    bval.mcmc <- proc.time()
    model.mcmc <- try(MCMCglmm(cbind(vy1, vy2) ~ trait - 1 + 
                                 at.level(trait, 1):(mX11+mX12+mX13+mX14)+ 
                                 at.level(trait, 2):(mX21+mX22), 
                               random = ~us(trait):factor, 
                               rcov = ~idh(trait):units, 
                               family = c("poisson", "poisson"), 
                               data = df.mcmc, 
                               prior = prior, 
                               nitt = 13000, 
                               burnin = 3000, 
                               thin = 10, 
                               verbose = FALSE), 
                      silent = TRUE)
    eval.mcmc <- proc.time()
    
    if((class(model.mcmc)[1])!="try-error"){ 
      converged.mcmc <- TRUE
      time.mcmc <- eval.mcmc[3] - bval.mcmc[3]
      beta.mcmc[[1]] <- as.numeric(summary(model.mcmc)$solutions[c(1,3:6),1])
      beta.mcmc[[2]] <- as.numeric(summary(model.mcmc)$solutions[c(2,7:8),1])
      sigma.mcmc[1] <- as.numeric(summary(model.mcmc)$Gcovariances[1,1])
      sigma.mcmc[2] <- as.numeric(summary(model.mcmc)$Gcovariances[4,1])
      rho.mcmc <- as.numeric(summary(model.mcmc)$Gcovariances[2,1])/sqrt(sigma.mcmc[1]*sigma.mcmc[2])
    } else {
      converged.mcmc <- FALSE
    }
    
    res.mcmc[[ITER,trial]] <- list(converged=converged.mcmc, time=time.mcmc, beta=beta.mcmc, sigma=sigma.mcmc, rho=rho.mcmc)
    
    
  }
}

################ result #######

index0.gvam <- list()
index0.mcmc <- list()
for (ITER in 1:length(M)) {
  index0.gvam[[ITER]] <- which(sapply(res.gvam[ITER,], function(list) all(list$converged)&!is.nan(list$rho)))
  index0.mcmc[[ITER]] <- which(sapply(res.mcmc[ITER,], function(list) list$converged))
}

meantime.gvam <- numeric(length(M))
meantime.mcmc <- numeric(length(M))
for (ITER in 1:length(M)) {
  if (length(index0.gvam[[ITER]]) > 0) {
    meantime.gvam[ITER] <- mean(sapply(res.gvam[ITER, index0.gvam[[ITER]]], function(list) list$time))
  }
  if (length(index0.mcmc[[ITER]]) > 0) {
    meantime.mcmc[ITER] <- mean(sapply(res.mcmc[ITER, index0.mcmc[[ITER]]], function(list) list$time))
  }
}

beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) list()))
beta.mcmc <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) list()))
mean.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
mean.beta.mcmc <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
bias.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
bias.beta.mcmc <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
mse.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
mse.beta.mcmc <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
for (k in 1:K) {
  for (p in 1:(P[k]+1)) {
    for (ITER in 1:length(M)) {
      if (length(index0.gvam[[ITER]]) > 0) {
        beta.gvam[[k]][[p]][[ITER]] <- sapply(res.gvam[ITER, index0.gvam[[ITER]]], function(list) list$beta[[k]][p])
        mean.beta.gvam[[k]][[p]][ITER] <- mean(beta.gvam[[k]][[p]][[ITER]])
        bias.beta.gvam[[k]][[p]][ITER] <- (mean.beta.gvam[[k]][[p]][ITER]-true.beta[[k]][p])/abs(true.beta[[k]][p])
        mse.beta.gvam[[k]][[p]][ITER] <- mean((beta.gvam[[k]][[p]][[ITER]]-true.beta[[k]][p])^2)
      }
      if (length(index0.mcmc[[ITER]]) > 0) {
        beta.mcmc[[k]][[p]][[ITER]] <- sapply(res.mcmc[ITER, index0.mcmc[[ITER]]], function(list) list$beta[[k]][p])
        mean.beta.mcmc[[k]][[p]][ITER] <- mean(beta.mcmc[[k]][[p]][[ITER]])
        bias.beta.mcmc[[k]][[p]][ITER] <- (mean.beta.mcmc[[k]][[p]][ITER]-true.beta[[k]][p])/abs(true.beta[[k]][p])
        mse.beta.mcmc[[k]][[p]][ITER] <- mean((beta.mcmc[[k]][[p]][[ITER]]-true.beta[[k]][p])^2)
      }
    }
  }
}      

sigma.gvam <- lapply(1:K, function(k) list())
sigma.mcmc <- lapply(1:K, function(k) list())
mean.sigma.gvam <- lapply(1:K, function(k)  numeric(length(M)))
mean.sigma.mcmc <- lapply(1:K, function(k)  numeric(length(M)))
bias.sigma.gvam <- lapply(1:K, function(k)  numeric(length(M)))
bias.sigma.mcmc <- lapply(1:K, function(k)  numeric(length(M)))
mse.sigma.gvam <- lapply(1:K, function(k)  numeric(length(M)))
mse.sigma.mcmc <- lapply(1:K, function(k)  numeric(length(M)))
for (k in 1:K) {
  for (ITER in 1:length(M)) {
    if (length(index0.gvam[[ITER]]) > 0) {
      sigma.gvam[[k]][[ITER]] <- sapply(res.gvam[ITER, index0.gvam[[ITER]]], function(list) list$sigma[k])
      mean.sigma.gvam[[k]][ITER] <- mean(sigma.gvam[[k]][[ITER]])
      bias.sigma.gvam[[k]][ITER] <- (mean.sigma.gvam[[k]][ITER]-true.Sigma[k,k])/abs(true.Sigma[k,k])
      mse.sigma.gvam[[k]][ITER] <- mean((sigma.gvam[[k]][[ITER]]-true.Sigma[k,k])^2)
    }
    if (length(index0.mcmc[[ITER]]) > 0) {
      sigma.mcmc[[k]][[ITER]] <- sapply(res.mcmc[ITER, index0.mcmc[[ITER]]], function(list) list$sigma[k])
      mean.sigma.mcmc[[k]][ITER] <- mean(sigma.mcmc[[k]][[ITER]])
      bias.sigma.mcmc[[k]][ITER] <- (mean.sigma.mcmc[[k]][ITER]-true.Sigma[k,k])/abs(true.Sigma[k,k])
      mse.sigma.mcmc[[k]][ITER] <- mean((sigma.mcmc[[k]][[ITER]]-true.Sigma[k,k])^2)
    }
  }
}

rho.gvam <- list()
rho.mcmc <- list()
mean.rho.gvam <- numeric(length(M))
mean.rho.mcmc <- numeric(length(M))
bias.rho.gvam <- numeric(length(M))
bias.rho.mcmc <- numeric(length(M))
mse.rho.gvam <- numeric(length(M))
mse.rho.mcmc <- numeric(length(M))
for (ITER in 1:length(M)) {
  
  if (length(index0.gvam[[ITER]]) > 0) {
    rho.gvam[[ITER]] <- sapply(res.gvam[ITER, index0.gvam[[ITER]]], function(list) list$rho)
    mean.rho.gvam[ITER] <- mean(rho.gvam[[ITER]])
    bias.rho.gvam[ITER] <- (mean.rho.gvam[ITER]-true.Rho[1,2])/abs(true.Rho[1,2])
    mse.rho.gvam[ITER] <- mean((rho.gvam[[ITER]]-true.Rho[1,2])^2)
  }
  
  if (length(index0.mcmc[[ITER]]) > 0) {
    rho.mcmc[[ITER]] <- sapply(res.mcmc[ITER, index0.mcmc[[ITER]]], function(list) list$rho)
    mean.rho.mcmc[ITER] <- mean(rho.mcmc[[ITER]])
    bias.rho.mcmc[ITER] <- (mean.rho.mcmc[ITER]-true.Rho[1,2])/abs(true.Rho[1,2])
    mse.rho.mcmc[ITER] <- mean((rho.mcmc[[ITER]]-true.Rho[1,2])^2)
  }
}


df.time <- round(as.data.frame(rbind(meantime.gvam, meantime.mcmc)), 2)
colnames(df.time) <- paste0('m=',1:length(M))
df.estimator <- list()
for (ITER in 1:length(M)){
  df.estimator[[paste0('m=',M[ITER])]] <- round(data.frame(
    true.value = round(c(unlist(true.beta), diag(true.Sigma), true.Rho[1,2]),2),
    bias.gvam = c(
      unlist(lapply(bias.beta.gvam, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(bias.sigma.gvam, function(l){l[ITER]})),
      bias.rho.gvam[ITER]),
    bias.mcmc = c(
      unlist(lapply(bias.beta.mcmc, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(bias.sigma.mcmc, function(l){l[ITER]})),
      bias.rho.mcmc[ITER]),
    mse.gvam = c(
      unlist(lapply(mse.beta.gvam, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(mse.sigma.gvam, function(l){l[ITER]})),
      mse.rho.gvam[ITER]),
    mse.mcmc = c(
      unlist(lapply(mse.beta.mcmc, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(mse.sigma.mcmc, function(l){l[ITER]})),
      mse.rho.mcmc[ITER])
  ), 5)
  rownames(df.estimator[[paste0('m=',M[ITER])]]) <- c(
    paste0('beta1_',0:(length(true.beta[[1]])-1)),
    paste0('beta2_',0:(length(true.beta[[2]])-1)),
    'sigma1^2', 'sigma2^2', 'rho')
}
print(df.time)
print(df.estimator)



