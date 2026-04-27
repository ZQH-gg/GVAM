rm(list = ls())
getwd()
setwd("C:/Users/admin/Desktop/gvam_code")
library(scales)
library(ggplot2)
library(Matrix)
library(mvtnorm) 
library(magic)
source("CalculateB.Rs")
source("GVArandint.Rs")
source("mglmmVB.R")
################ functions ##############
## Dataset generation for independent cluster GLMM; u_i ~ N(0,Sigma)
## Intercept must be manually included in fixed effects and random effects if desigray
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

GVAM.FIT <- function(Xdependent, id, Y, fixed.MM, beta0, gamma0, mu0, zeta0, family,MAXITER,EPS_GRAD,EPS_PAR) {
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
  ############# gva ###############
  unique.ids <- unique(id)
  m  <- length(unique.ids)
  gamma <- gamma0 
  vmu   <- matrix(mu0,m,1)
  vzeta <- matrix(zeta0,m,1)
  for (k in 1:K){
    vy <- Y[,k]
    mX <- fixed.MM[[k]]	
    MFVBeta <- matrix(beta0,p[[k]],1)
    
    cat("Fitting model using GVA", "  k=", k, "\n")						
    bval  <- proc.time()
    res <- try(GVA.fit <- GVA.randint.FIT(id,vy,mX,MFVBeta,sigma2=exp(gamma),vmu,vlambda=exp(vzeta),
                                          family,MAXITER,EPS_GRAD,EPS_PAR), silent = TRUE) 
    eval  <- proc.time()
    
    if((class(res)[1])!="try-error" && res$iterations<MAXITER){ 
      converged[k] <- TRUE
      times[k] <- eval[3] - bval[3]
      beta[[k]] <- res$vbeta
      sigma[k] <- res$sigma2
    } else {
      converged[k] <- FALSE
    }
  }
  
  ############# moment #############
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
    times[K+1] <- eval[3] - bval[3]
    names(rho) <- sapply(combinations, function(x) paste0("rho", x[1], x[2]))
  }
  out <- list(converged=converged, times=times, beta=beta, sigma=sigma, rho=rho)
  return(out)
}

################ simulation setting #####################
Xdependent=FALSE 
M <- c(10,50,100,500,1000,3000,5000,9000) ## Number of subjects  
ID <- lapply(M, function(m) rep(1:m, round(rnorm(m, mean = 10, sd = 1))))
MS <- 200 ## number of simulations
K = 2 ## number of responses
P = c(4, 2) ## numbers of fixed effects including intercept for different responses
## beta including intercept
true.beta <- list(
  beta1 = c(-1.0, 0.5, -0.8, 0.2, 0.6), 
  beta2 = c(-0.5, 0.3, -1.0)  
)  ## include intercept
true.Sigma <- matrix(c(1,0.5,0.5,0.8),ncol=2)
true.Rho <- true.Sigma/sqrt(diag(true.Sigma) %o% diag(true.Sigma))
combinations <- combn(K, 2, simplify = FALSE)

################ begin simulating #####################
res.gvam <- array(list(), dim = c(length(M), MS))
res.MFVB <- array(list(), dim = c(length(M), MS))
ITER0 <- 1
for (ITER in ITER0:length(M)) {
  m <- M[ITER]
  id <- ID[[ITER]] 
  trial0 <- 1
  for (trial in trial0:MS) {
    set.seed((ITER-1)*MS + trial)
    ########## data preparation(id, Y, fixed.MM) ########
 
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
    
    ########## gvam fitting #########
    Xdependent <- Xdependent
    beta0 <- 0  ## fixed effect
    gamma0 <- 0  ## log(variance) of random effect
    mu0   <- 0  ## expectaion of variational distribution
    zeta0 <- 0  ## log(variance) of variational distribution
    family = "POISSON"
    MAXITER = 100
    EPS_GRAD = 1e-5
    EPS_PAR = 1e-5
    cat("m=",m,"    trial=",trial, "\n")
    res.gvam[[ITER,trial]] <- GVAM.FIT(Xdependent=Xdependent, id=id, Y=Y, fixed.MM=fixed.MM, beta0=beta0, gamma0=gamma0, mu0=mu0, zeta0=zeta0, family = family,MAXITER,EPS_GRAD,EPS_PAR)
    
    ########## MFVB fitting #########
    converged.MFVB <- FALSE
    time.MFVB <- 0
    beta.MFVB <- lapply(1:K, function(k) rep(0,P[k]))
    sigma.MFVB <- rep(0,K)
    rho.MFVB <- rep(0,length(combinations))
    
    y_MFVB <- as.data.frame(Y)
    x_MFVB <- fixed.MM   
    z_MFVB <- ran.MM        
    id_vec <- id              
    responseType <- rep("Poisson", K)
    stopifnot(
      nrow(y_MFVB) == length(id_vec),
      nrow(y_MFVB) == nrow(x_MFVB[[1]]),
      nrow(y_MFVB) == nrow(z_MFVB[[1]]),
      length(responseType) == ncol(y_MFVB)
    )
    MAXITER = 100
    
    t.MFVB <- system.time({
      model.MFVB <- mglmmVB(
        y = y_MFVB,
        x = x_MFVB,
        z = z_MFVB,
        id = id_vec,
        responseType = responseType,
        method = "KMW",          
        doStreamlined = TRUE,
        maxIter = MAXITER,           
        useMatForDg = TRUE
      )
    })
    
    if((class(model.MFVB)[1])!="try-error" && model.MFVB$itnum<MAXITER){ 
      converged.MFVB <- TRUE
      time.MFVB <- t.MFVB[3]
      
      MFVBeta.MFVB <- model.MFVB$mu.q.betauG
      pr <- sapply(fixed.MM, ncol)   
      cum_pr <- c(0, cumsum(pr))
      beta.MFVB <- lapply(1:K, function(k) MFVBeta.MFVB[(cum_pr[k]+1):cum_pr[k+1]])
      
      Sigma.MFVB <- solve(model.MFVB$M.q.inv.SigmaR)
      sigma.MFVB <- diag(Sigma.MFVB)
      Rho.MFVB <- Sigma.MFVB/sqrt(sigma.MFVB%o%sigma.MFVB)
      rho.MFVB <- Rho.MFVB[upper.tri(Rho.MFVB, diag = FALSE)]
      names(rho.MFVB) <- sapply(combinations, function(x) paste0("rho", x[1], x[2]))
    } else {
      converged.MFVB <- FALSE
    }
    
    res.MFVB[[ITER,trial]] <- list(converged=converged.MFVB, time=time.MFVB, beta=beta.MFVB, sigma=sigma.MFVB, rho=rho.MFVB)
    
    ########## result saving #########
    save.image("C:/Users/admin/Desktop/gvam_code/simulation.RData")
  }
}

################ result analysis #######
#load("C:/Users/admin/Desktop/gvam_code/simulation.RData")
res.MFVB<-res.VB
for (ITER in 1:length(M)) {
  for (trial in 1:MS) {
    sigma<-res.MFVB[[ITER,trial]]$sigma
    rho<-res.MFVB[[ITER,trial]]$rho/sqrt(sigma[1]*sigma[2])
    res.MFVB[[ITER,trial]]$rho <- rho
  }
}
# converged index ----------
is_outlier <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  x < (q[1] - 1.5 * iqr) | x > (q[2] + 1.5 * iqr)
}

common_index0 <- list()  
final_index    <- list() 
for (ITER in 1:length(M)) {
  idx_g <- which(sapply(res.gvam[ITER,], function(l) all(l$converged) & all(l$rho>-1&l$rho<1)))
  idx_m <- which(sapply(res.MFVB[ITER,], function(l) l$converged & all(l$rho>-1&l$rho<1)))
  common_index0[[ITER]] <- intersect(idx_g, idx_m)
  
  idx <- common_index0[[ITER]]
  n <- length(idx)
  if (n == 0) {
    final_index[[ITER]] <- integer(0)
    next
  }
  
  params_g <- lapply(res.gvam[ITER, idx], function(o) c(unlist(o$beta), o$sigma, unlist(o$rho)))
  mat_g    <- do.call(rbind, params_g)
  out_g    <- apply(mat_g, 2, is_outlier)
  out_g    <- rowSums(out_g) > 0  
  
  params_m <- lapply(res.MFVB[ITER, idx], function(o) c(unlist(o$beta), o$sigma, unlist(o$rho)))
  mat_m    <- do.call(rbind, params_m)
  out_m    <- apply(mat_m, 2, is_outlier)
  out_m    <- rowSums(out_m) > 0  
  
  keep <- !out_g & !out_m
  final_index[[ITER]] <- idx[keep]
}
index <- lapply(final_index,function(x) x[1:100])
print(sapply(index,length))
# statistics -------
# Computational Cost
meantime.gvam <- numeric(length(M))
meantime.MFVB <- numeric(length(M))
for (ITER in 1:length(M)) {
  if (length(index[[ITER]]) > 0) {
    meantime.gvam[ITER] <- mean(sapply(res.gvam[ITER, index[[ITER]]], function(list) sum(list$times)))
    meantime.MFVB[ITER] <- mean(sapply(res.MFVB[ITER, index[[ITER]]], function(list) list$time))
  }
}
# beta.hat
beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) list()))
beta.MFVB <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) list()))
mean.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
mean.beta.MFVB <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
relative_bias.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
relative_bias.beta.MFVB <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
mse.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
mse.beta.MFVB <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(length(M))))
for (k in 1:K) {
  for (p in 1:(P[k]+1)) {
    for (ITER in 1:length(M)) {
      if (length(index[[ITER]]) > 0) {
        beta.gvam[[k]][[p]][[ITER]] <- sapply(res.gvam[ITER, index[[ITER]]], function(list) list$beta[[k]][p])
        mean.beta.gvam[[k]][[p]][ITER] <- mean(beta.gvam[[k]][[p]][[ITER]])
        relative_bias.beta.gvam[[k]][[p]][ITER] <- (mean.beta.gvam[[k]][[p]][ITER]-true.beta[[k]][p])/abs(true.beta[[k]][p])
        mse.beta.gvam[[k]][[p]][ITER] <- mean((beta.gvam[[k]][[p]][[ITER]]-true.beta[[k]][p])^2)
        
        beta.MFVB[[k]][[p]][[ITER]] <- sapply(res.MFVB[ITER, index[[ITER]]], function(list) list$beta[[k]][p])
        mean.beta.MFVB[[k]][[p]][ITER] <- mean(beta.MFVB[[k]][[p]][[ITER]])
        relative_bias.beta.MFVB[[k]][[p]][ITER] <- (mean.beta.MFVB[[k]][[p]][ITER]-true.beta[[k]][p])/abs(true.beta[[k]][p])
        mse.beta.MFVB[[k]][[p]][ITER] <- mean((beta.MFVB[[k]][[p]][[ITER]]-true.beta[[k]][p])^2)
      }
    }
  }
}      
# sigma.hat
sigma.gvam <- lapply(1:K, function(k) list())
sigma.MFVB <- lapply(1:K, function(k) list())
mean.sigma.gvam <- lapply(1:K, function(k)  numeric(length(M)))
mean.sigma.MFVB <- lapply(1:K, function(k)  numeric(length(M)))
relative_bias.sigma.gvam <- lapply(1:K, function(k)  numeric(length(M)))
relative_bias.sigma.MFVB <- lapply(1:K, function(k)  numeric(length(M)))
mse.sigma.gvam <- lapply(1:K, function(k)  numeric(length(M)))
mse.sigma.MFVB <- lapply(1:K, function(k)  numeric(length(M)))
for (k in 1:K) {
  for (ITER in 1:length(M)) {
    if (length(index[[ITER]]) > 0) {
      sigma.gvam[[k]][[ITER]] <- sapply(res.gvam[ITER, index[[ITER]]], function(list) list$sigma[k])
      mean.sigma.gvam[[k]][ITER] <- mean(sigma.gvam[[k]][[ITER]])
      relative_bias.sigma.gvam[[k]][ITER] <- (mean.sigma.gvam[[k]][ITER]-true.Sigma[k,k])/abs(true.Sigma[k,k])
      mse.sigma.gvam[[k]][ITER] <- mean((sigma.gvam[[k]][[ITER]]-true.Sigma[k,k])^2)
      
      sigma.MFVB[[k]][[ITER]] <- sapply(res.MFVB[ITER, index[[ITER]]], function(list) list$sigma[k])
      mean.sigma.MFVB[[k]][ITER] <- mean(sigma.MFVB[[k]][[ITER]])
      relative_bias.sigma.MFVB[[k]][ITER] <- (mean.sigma.MFVB[[k]][ITER]-true.Sigma[k,k])/abs(true.Sigma[k,k])
      mse.sigma.MFVB[[k]][ITER] <- mean((sigma.MFVB[[k]][[ITER]]-true.Sigma[k,k])^2)
    }
  }
}
# rho.hat
rho.gvam <- list()
rho.MFVB <- list()
mean.rho.gvam <- list()
mean.rho.MFVB <- list()
relative_bias.rho.gvam <- list()
relative_bias.rho.MFVB <- list()
mse.rho.gvam <- list()
mse.rho.MFVB <- list()
for (x in combinations) {
  rho_name<-paste0("rho", x[1], x[2])
  for (ITER in 1:length(M)) {
    if (length(index[[ITER]]) > 0) {
      rho.gvam[[rho_name]][[ITER]] <- sapply(res.gvam[ITER, index[[ITER]]], function(list) list$rho[[rho_name]])
      mean.rho.gvam[[rho_name]][ITER] <- mean(rho.gvam[[rho_name]][[ITER]])
      relative_bias.rho.gvam[[rho_name]][ITER] <- (mean.rho.gvam[[rho_name]][ITER]-true.Rho[x[1],x[2]])/abs(true.Rho[x[1],x[2]])
      mse.rho.gvam[[rho_name]][ITER] <- mean((rho.gvam[[rho_name]][[ITER]]-true.Rho[x[1],x[2]])^2)
      
      rho.MFVB[[rho_name]][[ITER]] <- sapply(res.MFVB[ITER, index[[ITER]]], function(list) list$rho[[rho_name]])
      mean.rho.MFVB[[rho_name]][ITER] <- mean(rho.MFVB[[rho_name]][[ITER]])
      relative_bias.rho.MFVB[[rho_name]][ITER] <- (mean.rho.MFVB[[rho_name]][ITER]-true.Rho[x[1],x[2]])/abs(true.Rho[x[1],x[2]])
      mse.rho.MFVB[[rho_name]][ITER] <- mean((rho.MFVB[[rho_name]][[ITER]]-true.Rho[x[1],x[2]])^2)
    }
  }
}


# plot of Computational Cost -----------
data_meantime.gvam <- data.frame(x = sapply(ID,length), y = meantime.gvam)
data_meantime.gvam$log_x <- log10(data_meantime.gvam$x)
data_meantime.gvam$log_y <- log10(data_meantime.gvam$y)
data_meantime.gvam$method <- "GVA-M"  

data_meantime.MFVB <- data.frame(x = sapply(ID,length), y = meantime.MFVB)
data_meantime.MFVB$log_x <- log10(data_meantime.MFVB$x)
data_meantime.MFVB$log_y <- log10(data_meantime.MFVB$y)
data_meantime.MFVB$method <- "MFVB"  


data_combined <- rbind(data_meantime.gvam, data_meantime.MFVB)

lm_model.gvam <- lm(log_y ~ log_x, data = data_meantime.gvam)
slope.gvam <- summary(lm_model.gvam)$coefficients[2,1]
lm_model.MFVB <- lm(log_y ~ log_x, data = data_meantime.MFVB)
slope.MFVB <- summary(lm_model.MFVB)$coefficients[2,1]

plot_time <- ggplot(data_combined, aes(x = log_x, y = log_y, color = method)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  
  annotate("text", 
           x = max(data_meantime.gvam$log_x) * 0.8, 
           y = max(data_meantime.gvam$log_y) * 0.8, 
           label = paste("Slope:", round(slope.gvam, 2)), 
           color = "blue", size = 4) +
  annotate("text", 
           x = max(data_meantime.MFVB$log_x) * 0.8, 
           y = max(data_meantime.MFVB$log_y) * 0.8, 
           label = paste("Slope:", round(slope.MFVB, 2)), 
           color = "darkgreen", size = 4) +
  
  scale_color_manual(values = c("GVA-M" = "blue", "MFVB" = "darkgreen")) +
  labs(x = "Log (Sample Size)", y = "Log (Time/s)", title = "CPU sec.", color = "Method") +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0, size = 14),  
    axis.title.x = element_text(size = 12),          
    axis.text.x = element_text(size = 10),            
    axis.title.y = element_text(size = 12),           
    axis.text.y = element_text(size = 10),             
    legend.title = element_text(size = 12),     
    legend.text = element_text(size = 10)      
  )

print(plot_time)

# box-plots ------------
ylim_beta <- list(
  list(c(-2.14,0.03),c(0.15,0.83),c(-1.2,-0.5),c(-0.15,0.6),c(0.3,9.5)),
  list(c(-1.4,0.2),c(0,0.6),c(-1.3,-0.7))
)
ylim_sigma <- list(
  c(0,3.5),
  c(0,3.2)
)
# boxplot_beta-------
boxplot_beta.gvam <- lapply(1:K, function(k) list())
boxplot_beta.MFVB <- lapply(1:K, function(k) list())
for (k in 1:K) {
  for (p in 1:(P[k]+1)) {
    # ---------- GVA-M ----------
    data_beta.gvam <- data.frame(
      group = factor(rep(M, sapply(index, length)), levels = M),
      value = unlist(lapply(1:length(M), function(ITER) beta.gvam[[k]][[p]][[ITER]]))
    )
    boxplot_beta.gvam[[k]][[p]] <- ggplot(data_beta.gvam, aes(x = group, y = value, fill = group)) +
      geom_boxplot(fill = "gray70", color = "gray30", outlier.shape = NA) +
      geom_hline(yintercept = true.beta[[k]][p], color = "red", linetype = "solid", size = 0.8) +
      labs(x = "Number of Subjects", y = NULL, title = substitute(beta[k*p]*" GVA-M", list(k = k, p = p - 1))) +
      ylim(ylim_beta[[k]][[p]][1],ylim_beta[[k]][[p]][2]) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0, size = 16),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )
    
    # ---------- MFVB ----------
    data_beta.MFVB <- data.frame(
      group = factor(rep(M, sapply(index, length)), levels = M),
      value = unlist(lapply(1:length(M), function(ITER) beta.MFVB[[k]][[p]][[ITER]]))
    )
    boxplot_beta.MFVB[[k]][[p]] <- ggplot(data_beta.MFVB, aes(x = group, y = value, fill = group)) +
      geom_boxplot(fill = "gray70", color = "gray30", outlier.shape = NA) +
      geom_hline(yintercept = true.beta[[k]][p], color = "red", linetype = "solid", size = 0.8) +
      labs(x = "Number of Subjects", y = NULL, title = substitute(beta[k*p]*" MFVB", list(k = k, p = p - 1))) +
      ylim(ylim_beta[[k]][[p]][1],ylim_beta[[k]][[p]][2]) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0, size = 16),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )
    
  }
}
print(boxplot_beta.gvam)
print(boxplot_beta.MFVB)
# boxplot_sigma--------
boxplot_sigma.gvam <- list()
boxplot_sigma.MFVB <- list()
for (k in 1:K) {
  # ---------- GVA-M ----------
  data_sigma.gvam <- data.frame(
    group = factor(rep(M, sapply(index, length)), levels = M),
    value = unlist(lapply(1:length(M), function(ITER) sigma.gvam[[k]][[ITER]]))
  )
  boxplot_sigma.gvam[[k]] <- ggplot(data_sigma.gvam, aes(x = group, y = value, fill = group)) +
    geom_boxplot(fill = "gray70", color = "gray30", outlier.shape = NA) +
    geom_hline(yintercept = true.Sigma[k, k], color = "red", linetype = "solid", size = 0.8) +
    labs(x = "Number of Subjects", y = NULL, title = substitute(sigma[k]*" GVA-M", list(k = k))) +
    ylim(ylim_sigma[[k]][1],ylim_sigma[[k]][2]) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 16),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  # ---------- MFVB ----------
  data_sigma.MFVB <- data.frame(
    group = factor(rep(M, sapply(index, length)), levels = M),
    value = unlist(lapply(1:length(M), function(ITER) sigma.MFVB[[k]][[ITER]]))
  )
  boxplot_sigma.MFVB[[k]] <- ggplot(data_sigma.MFVB, aes(x = group, y = value, fill = group)) +
    geom_boxplot(fill = "gray70", color = "gray30", outlier.shape = NA) +
    geom_hline(yintercept = true.Sigma[k, k], color = "red", linetype = "solid", size = 0.8) +
    labs(x = "Number of Subjects", y = NULL, title = substitute(sigma[k]*" MFVB", list(k = k))) +
    ylim(ylim_sigma[[k]][1],ylim_sigma[[k]][2]) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 16),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}
print(boxplot_sigma.gvam)
print(boxplot_sigma.MFVB)
# boxplot_rho ----------
boxplot_rho.gvam <- list()
boxplot_rho.MFVB <- list()
for (x in combinations) {
  rho_name <- paste0("rho", x[1], x[2])
  # ---------- GVA-M ----------
  data_rho.gvam <- data.frame(
    group = factor(rep(M, sapply(index, length)), levels = M),
    value = unlist(lapply(1:length(M), function(ITER) rho.gvam[[rho_name]][[ITER]]))
  )
  boxplot_rho.gvam[[rho_name]] <- ggplot(data_rho.gvam, aes(x = group, y = value, fill = group)) +
    geom_boxplot(fill = "gray70", color = "gray30", outlier.shape = NA) +
    geom_hline(yintercept = true.Rho[x[1], x[2]], color = "red", linetype = "solid", size = 0.8) +
    labs(x = "Number of Subjects", y = NULL, title = substitute(rho[x]*" GVA-M", list(x = paste0(x[1], x[2])))) +
    ylim(-0.9,1)+
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 16),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  # ---------- MFVB ----------
  data_rho.MFVB <- data.frame(
    group = factor(rep(M, sapply(index, length)), levels = M),
    value = unlist(lapply(1:length(M), function(ITER) rho.MFVB[[rho_name]][[ITER]]))
  )
  boxplot_rho.MFVB[[rho_name]] <- ggplot(data_rho.MFVB, aes(x = group, y = value, fill = group)) +
    geom_boxplot(fill = "gray70", color = "gray30", outlier.shape = NA) +
    geom_hline(yintercept = true.Rho[x[1], x[2]], color = "red", linetype = "solid", size = 0.8) +
    labs(x = "Number of Subjects", y = NULL, title = substitute(rho[x]*" MFVB", list(x = paste0(x[1], x[2])))) +
    ylim(-0.9,1)+
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 16),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}
print(boxplot_rho.gvam)
print(boxplot_rho.MFVB)
# display statistics -----------
df.time <- round(as.data.frame(rbind(meantime.gvam, meantime.MFVB)), 2)
colnames(df.time) <- paste0('m=',1:length(M))
df.estimator <- list()
for (ITER in 1:length(M)){
  df.estimator[[paste0('m=',M[ITER])]] <- round(data.frame(
    true.value = round(c(unlist(true.beta), diag(true.Sigma), true.Rho[upper.tri(true.Rho, diag = FALSE)]),2),
    relative_bias.gvam = c(
      unlist(lapply(relative_bias.beta.gvam, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(relative_bias.sigma.gvam, function(l){l[ITER]})),
      unlist(lapply(relative_bias.rho.gvam, function(l){l[ITER]}))
    ),
    mse.gvam = c(
      unlist(lapply(mse.beta.gvam, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(mse.sigma.gvam, function(l){l[ITER]})),
      unlist(lapply(mse.rho.gvam, function(l){l[ITER]}))
    ),
    relative_bias.MFVB = c(
      unlist(lapply(relative_bias.beta.MFVB, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(relative_bias.sigma.MFVB, function(l){l[ITER]})),
      unlist(lapply(relative_bias.rho.MFVB, function(l){l[ITER]}))
    ),
    mse.MFVB = c(
      unlist(lapply(mse.beta.MFVB, function(list){
        sapply(list, function(l){l[ITER]})})),
      unlist(lapply(mse.sigma.MFVB, function(l){l[ITER]})),
      unlist(lapply(mse.rho.MFVB, function(l){l[ITER]}))
    )
  ), 5)
  rownames(df.estimator[[paste0('m=',M[ITER])]]) <- c(
    unlist(lapply(1:K, function(k) paste0('beta',k,'_',0:(length(true.beta[[k]])-1)))),
    paste0('sigma',1:K,'^2'), 
    sapply(combinations, function(x) paste0("rho", x[1], x[2]))
  )
}
print(df.time)
print(df.estimator)



###############



