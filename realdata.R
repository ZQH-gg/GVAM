rm(list = ls())
getwd()
setwd("C:/Users/admin/Desktop/gvam_code")
library(psych)
library(dplyr)
library(Matrix)
library(MCMCglmm)
library(MASS)
library(ggplot2)
library(haven)
library(glmmTMB)
source("CalculateB.Rs")
source("GVArandint.Rs")
########### functions ##########
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

########### real data ################
df_bac <- read.csv("taxonomy_counts_g.csv")
df_nut <- read.table("nutrition_totals.txt", header = TRUE, sep = "")
set.seed(8)
rownames(df_bac) <- df_bac[,1]
df_bac <- df_bac[,-1]
df_bac <- as.data.frame(t(df_bac))

df1 <- data.frame(
  X.SampleID = rownames(df_bac),
  y1 = df_bac[, grep("Clostridioides", colnames(df_bac), ignore.case = TRUE)],
  y2 = df_bac[, grep("Enterococcus", colnames(df_bac), ignore.case = TRUE)]
)
df2 <- df_nut[,c(1,2,5:69)]
df_merge <- merge(df1, df2, by = "X.SampleID")
attr(df_merge$y1, "label") <- "Clostridioides"
attr(df_merge$y2, "label") <- "Enterococcus"

data<-df_merge
data$y1[data$y1==0]<-1
data$y2[data$y2==0]<-1

y1_log<-log(data$y1)
full_model1 <- lm(y1_log ~ ., data = cbind(y1_log,data[,-(1:4)]))
stepwise_model1 <- step(full_model1, direction = "both",k = qchisq(0.001, 1, lower.tail = FALSE))
vars1 <-rownames(summary(stepwise_model1)$coefficients)[-1]

y2_log<-log(data$y2)
full_model2 <- lm(y2_log ~ ., data = cbind(y2_log,data[,-(1:4)]))
stepwise_model2 <- step(full_model2, direction = "both",k = qchisq(0.001, 1, lower.tail = FALSE))
vars2 <-rownames(summary(stepwise_model2)$coefficients)[-1]

# MCMC-select 
id <- as.vector(data$UserName)
try_df.mcmc <- data.frame(
  vy1 = data$y1,
  mX11 = data$CHOLE, 
  mX12 = data$S040,
  mX13 = data$S080,
  mX14 = data$S120,
  mX15 = data$CHOLN, 
  mX16 = data$VITE_ADD,
  vy2 = data$y2,
  mX21 = data$NIAC,
  mX22 = data$P204,
  mX23 = data$B12_ADD, 
  factor = factor(id),
  ID = c(1:length(id))
)
prior = list(
  R = list(V = diag(2), nu = 0.002), 
  G = list(G1 = list(V = diag(2), nu = 0.002))
)
try_model.mcmc <- try(MCMCglmm(cbind(vy1, vy2) ~ trait - 1 + 
                                 at.level(trait, 1):(mX11+mX12+mX13+mX14+mX15+mX16)+ 
                                 at.level(trait, 2):(mX21+mX22+mX23), 
                               random = ~us(trait):factor, 
                               rcov = ~idh(trait):units, 
                               family = c("poisson", "poisson"), 
                               data = try_df.mcmc, 
                               prior = prior, 
                               nitt = 13000, 
                               burnin = 3000, 
                               thin = 10, 
                               verbose = FALSE), 
                      silent = TRUE)

summary(try_model.mcmc)

# data final 
id <- as.vector(data$UserName)
Y <- as.matrix(data[c("y1","y2")])
fixed.MM <- list(
  as.matrix(cbind(1, data[c("S040","S080","S120","VITE_ADD")])),
  as.matrix(cbind(1, data[c('NIAC','P204')]))
)
ran.MM <- replicate(2, matrix(rep(1,length(id)),ncol=1), simplify = FALSE)

# GVAM fitting 
Xdependent <- FALSE  
beta0 <- 0  
gamma0 <- 0  
mu0   <- 0  
zeta0 <- 0  
family = "POISSON"
MAXITER = 100
EPS_GRAD = 1e-5
EPS_PAR = 1e-5
res.gvam <- GVAM.FIT(Xdependent=Xdependent, id=id, Y=Y, fixed.MM=fixed.MM, beta0=beta0, gamma0=gamma0, mu0=mu0, zeta0=zeta0, family = family,MAXITER,EPS_GRAD,EPS_PAR)

# MCMC fiting 
converged.mcmc <- FALSE
time.mcmc <- 0
beta.mcmc <- list(rep(0,5), rep(0,3))
sigma.mcmc <- rep(0,2)
rho.mcmc <- 0

df.mcmc <- data.frame(
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

res.mcmc <- list(converged=converged.mcmc, time=time.mcmc, beta=beta.mcmc, sigma=sigma.mcmc, rho=rho.mcmc)
summary(model.mcmc)
print(res.gvam)
print(res.mcmc)
############### bootstrap  #########
unique_ids <- unique(id)
id_mapping <- lapply(unique_ids, function(x) which(id == x))
names(id_mapping) <- unique_ids
MS <- 200
res.gvam.resample <- array(list(), dim = MS)
res.mcmc.resample <- array(list(), dim = MS)
trial0 <- 1
for (trial in trial0:MS) {
  set.seed(trial)
  cat("trial=",trial, "\n")
  resampled_subjects <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
  
  resample_indices <- unlist(lapply(resampled_subjects, function(x) {
    id_mapping[[as.character(x)]]
  }))
  resample_ids <- rep(1:length(resampled_subjects), 
                      times = sapply(resampled_subjects, function(x) 
                        length(id_mapping[[as.character(x)]])))
  
  Y.resample <- Y[resample_indices, ]
  fixed.MM.resample <- list(
    fixed.MM[[1]][resample_indices, ],
    fixed.MM[[2]][resample_indices, ]
  )
  ############## GVAM fitting ##################
  Xdependent <- FALSE  
  beta0 <- 0  
  gamma0 <- 0 
  mu0   <- 0  
  zeta0 <- 0  
  family = "POISSON"
  MAXITER = 100
  EPS_GRAD = 1e-5
  EPS_PAR = 1e-5
  
  res.gvam.resample[[trial]] <- GVAM.FIT(Xdependent=Xdependent, id=resample_ids, Y=Y.resample, fixed.MM=fixed.MM.resample, beta0=beta0, gamma0=gamma0, mu0=mu0, zeta0=zeta0, family = family,MAXITER,EPS_GRAD,EPS_PAR)
  
  ############## MCMC fiting ##########
  converged.mcmc.resample <- FALSE
  time.mcmc.resample <- 0
  beta.mcmc.resample <- list(rep(0,5), rep(0,3))
  sigma.mcmc.resample <- rep(0,2)
  rho.mcmc.resample <- 0
  
  df.mcmc.resample <- data.frame(
    vy1 = Y.resample[,1],
    mX11 = fixed.MM.resample[[1]][,2],
    mX12 = fixed.MM.resample[[1]][,3],
    mX13 = fixed.MM.resample[[1]][,4],
    mX14 = fixed.MM.resample[[1]][,5],
    vy2 = Y.resample[,2],
    mX21 = fixed.MM.resample[[2]][,2],
    mX22 = fixed.MM.resample[[2]][,3],
    factor = factor(resample_ids),
    ID = c(1:length(resample_ids))
  )
  prior = list(
    R = list(V = diag(2), nu = 0.002), 
    G = list(G1 = list(V = diag(2), nu = 0.002))
  )
  bval.mcmc.resample <- proc.time()
  model.mcmc.resample <- try(MCMCglmm(cbind(vy1, vy2) ~ trait - 1 +
                                        at.level(trait, 1):(mX11+mX12+mX13+mX14)+
                                        at.level(trait, 2):(mX21+mX22),
                                      random = ~us(trait):factor,
                                      rcov = ~idh(trait):units,
                                      family = c("poisson", "poisson"),
                                      data = df.mcmc.resample,
                                      prior = prior,
                                      nitt = 13000,
                                      burnin = 3000,
                                      thin = 10,
                                      verbose = FALSE),
                             silent = TRUE)
  eval.mcmc.resample <- proc.time()
  
  if((class(model.mcmc.resample)[1])!="try-error"){
    converged.mcmc.resample <- TRUE
    time.mcmc.resample <- eval.mcmc.resample[3] - bval.mcmc.resample[3]
    beta.mcmc.resample[[1]] <- as.numeric(summary(model.mcmc.resample)$solutions[c(1,3:6),1])
    beta.mcmc.resample[[2]] <- as.numeric(summary(model.mcmc.resample)$solutions[c(2,7:8),1])
    sigma.mcmc.resample[1] <- as.numeric(summary(model.mcmc.resample)$Gcovariances[1,1])
    sigma.mcmc.resample[2] <- as.numeric(summary(model.mcmc.resample)$Gcovariances[4,1])
    rho.mcmc.resample <- as.numeric(summary(model.mcmc.resample)$Gcovariances[2,1])/sqrt(sigma.mcmc.resample[1]*sigma.mcmc.resample[2])
  } else {
    converged.mcmc <- FALSE
  }
  
  res.mcmc.resample[[trial]] <- list(converged=converged.mcmc.resample, time=time.mcmc.resample, beta=beta.mcmc.resample, sigma=sigma.mcmc.resample, rho=rho.mcmc.resample)
  
  #result saving --------------
  save.image("C:/Users/admin/Desktop/gvam_code/realdata.RData")
}
########### result analysis ############
#load("C:/Users/admin/Desktop/gvam_code/realdata.RData")
describe(df_merge$y1) ## Clostridioides
describe(df_merge$y2) ## Enterococcus
describe(df_merge[c("S040","S080","S120","VITE_ADD",'NIAC','P204')])
summary(model.mcmc)
print(res.gvam)
print(res.mcmc)
# converged index --------
is_outlier <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  x < (q[1] - 1.5 * iqr) | x > (q[2] + 1.5 * iqr)
}

common_index0 <- c()  
final_index    <- c() 

idx_g <- which(sapply(res.gvam.resample, function(l) all(l$converged) & all(l$rho>-1&l$rho<1)))
idx_m <- which(sapply(res.mcmc.resample, function(l) l$converged & all(l$rho>-1&l$rho<1)))
common_index0 <- intersect(idx_g, idx_m)

idx <- common_index0
n <- length(idx)
if (n == 0) {
  final_index <- integer(0)
  next
}

params_g <- lapply(res.gvam.resample[idx], function(o) c(unlist(o$beta), o$sigma, unlist(o$rho)))
mat_g    <- do.call(rbind, params_g)
out_g    <- apply(mat_g, 2, is_outlier)
out_g    <- rowSums(out_g) > 0  

params_m <- lapply(res.mcmc.resample[idx], function(o) c(unlist(o$beta), o$sigma, unlist(o$rho)))
mat_m    <- do.call(rbind, params_m)
out_m    <- apply(mat_m, 2, is_outlier)
out_m    <- rowSums(out_m) > 0  

keep <- !out_g & !out_m
final_index <- idx[keep]

index <- final_index[1:100]
print(length(common_index0))
print(length(final_index))
print(length(index))
# statistics -------
dp <- 5 
K <- ncol(Y)
P <- sapply(fixed.MM, ncol)-1
if (length(index) > 0) {
  # gvam
  # Computational Cost
  meantime.gvam <- round(mean(sapply(res.gvam.resample[index], function(list) list$time)), dp)
  # beta.hat
  beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(MS)))
  mean.beta.gvam <- lapply(1:K, function(k) numeric(P[k]+1))
  sd.beta.gvam <- lapply(1:K, function(k) numeric(P[k]+1))
  conf.beta.gvam <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(2))) 
  # sigma.hat
  sigma.gvam <- lapply(1:K, function(k) numeric(MS))
  mean.sigma.gvam <- rep(0,K)
  sd.sigma.gvam <- rep(0,K)
  conf.sigma.gvam <- lapply(1:K, function(k) numeric(2))
  for (k in 1:K) {
    for (p in 1:(P[k]+1)) {
      beta.gvam[[k]][[p]] <- sapply(res.gvam.resample[index], function(list) list$beta[[k]][p])
      mean.beta.gvam[[k]][[p]] <- round(mean(beta.gvam[[k]][[p]]), dp)
      sd.beta.gvam[[k]][[p]] <- round(sd(beta.gvam[[k]][[p]]), dp)
      conf.beta.gvam[[k]][[p]] <- round(quantile(beta.gvam[[k]][[p]], c(0.025, 0.975)), dp)
    }
    # sigma.hat
    sigma.gvam[[k]] <- sapply(res.gvam.resample[index], function(list) list$sigma[k])
    mean.sigma.gvam[k] <- round(mean(sigma.gvam[[k]]), dp)
    sd.sigma.gvam[k] <- round(sd(sigma.gvam[[k]]), dp)
    conf.sigma.gvam[[k]] <- round(quantile(sigma.gvam[[k]], c(0.025, 0.975)), dp)
  }   
  # rho.hat
  rho.gvam <- sapply(res.gvam.resample[index], function(list) list$rho)
  mean.rho.gvam <- round(mean(rho.gvam), dp)
  sd.rho.gvam <- round(sd(rho.gvam), dp)
  conf.rho.gvam <- round(quantile(rho.gvam, c(0.025, 0.975)), dp)
  
  # mcmc
  # Computational Cost
  meantime.mcmc <- round(mean(sapply(res.mcmc.resample[index], function(list) list$time)), dp)
  # beta.hat
  beta.mcmc <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(MS)))
  mean.beta.mcmc <- lapply(1:K, function(k) numeric(P[k]+1))
  sd.beta.mcmc <- lapply(1:K, function(k) numeric(P[k]+1))
  conf.beta.mcmc <- lapply(1:K, function(k) lapply(1:(P[k]+1), function(p) numeric(2))) 
  # sigma.hat
  sigma.mcmc <- lapply(1:K, function(k) numeric(MS))
  mean.sigma.mcmc <- rep(0,K)
  sd.sigma.mcmc <- rep(0,K)
  conf.sigma.mcmc <- lapply(1:K, function(k) numeric(2))
  for (k in 1:K) {
    for (p in 1:(P[k]+1)) {
      beta.mcmc[[k]][[p]] <- sapply(res.mcmc.resample[index], function(list) list$beta[[k]][p])
      mean.beta.mcmc[[k]][[p]] <- round(mean(beta.mcmc[[k]][[p]]), dp)
      sd.beta.mcmc[[k]][[p]] <- round(sd(beta.mcmc[[k]][[p]]), dp)
      conf.beta.mcmc[[k]][[p]] <- round(quantile(beta.mcmc[[k]][[p]], c(0.025, 0.975)), dp)
    }
    # sigma.hat
    sigma.mcmc[[k]] <- sapply(res.mcmc.resample[index], function(list) list$sigma[k])
    mean.sigma.mcmc[k] <- round(mean(sigma.mcmc[[k]]), dp)
    sd.sigma.mcmc[k] <- round(sd(sigma.mcmc[[k]]), dp)
    conf.sigma.mcmc[[k]] <- round(quantile(sigma.mcmc[[k]], c(0.025, 0.975)), dp)
  }   
  # rho.hat
  rho.mcmc <- sapply(res.mcmc.resample[index], function(list) list$rho)
  mean.rho.mcmc <- round(mean(rho.mcmc), dp)
  sd.rho.mcmc <- round(sd(rho.mcmc), dp)
  conf.rho.mcmc <- round(quantile(rho.mcmc, c(0.025, 0.975)), dp)
}
# box-plots ------------
# boxplot_beta
boxplot_beta <- lapply(1:K, function(k) list())
for (k in 1:K) {
  for (p in 1:(P[k]+1)) {
    data_beta <- data.frame(
      Group = factor(
        rep(c("GVA-M", "MCMC"), c(length(index),length(index))),
        levels = c("GVA-M", "MCMC")
      ),
      Estimates = c(beta.gvam[[k]][[p]], beta.mcmc[[k]][[p]])
    )
    
    boxplot_beta[[k]][[p]] <- ggplot(data_beta, aes(x = Group, y = Estimates, fill = Group)) +
      geom_boxplot(fill = "gray70", color = "gray30",  outlier.shape = NA) +
      labs(x = NULL,y = NULL,title = substitute(beta[k*p], list(k = k, p = p - 1))) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0, size = 16),  
        axis.title.x = element_text(size = 14),           
        axis.text.x = element_text(size = 12),            
        axis.text.y = element_text(size = 12)             
      )
  }
}
# boxplot_sigma
boxplot_sigma <- list()
for (k in 1:K) {
  data_sigma <- data.frame(
    Group = rep(c("GVA-M", "MCMC"), c(length(index),length(index))),
    Estimates = c(sigma.gvam[[k]], sigma.mcmc[[k]])
  )
  
  boxplot_sigma[[k]] <- ggplot(data_sigma, aes(x = Group, y = Estimates, fill = Group)) +
    geom_boxplot(fill = "gray70", color = "gray30",  outlier.shape = NA) +
    labs(x = NULL,y = NULL,title = substitute(sigma[k], list(k = k))) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 16),  
      axis.title.x = element_text(size = 14),           
      axis.text.x = element_text(size = 12),            
      axis.text.y = element_text(size = 12)             
    )
}
# boxplot_rho
data_rho <- data.frame(
  Group = rep(c("GVA-M", "MCMC"), c(length(index),length(index))),
  Estimates = c(rho.gvam, rho.mcmc)
)

boxplot_rho <- ggplot(data_rho, aes(x = Group, y = Estimates, fill = Group)) +
  geom_boxplot(fill = "gray70", color = "gray30",  outlier.shape = NA) +
  labs(x = NULL,y = NULL,title = expression(rho)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),  
    axis.title.x = element_text(size = 14),           
    axis.text.x = element_text(size = 12),            
    axis.text.y = element_text(size = 12)             
  )

print(boxplot_beta)
print(boxplot_sigma)
print(boxplot_rho)
# display statistics-----
df.time <- round(data.frame(
  rawtime.gvam=sum(res.gvam$times), 
  rawtime.mcmc=res.mcmc$time,
  meantime.gvam, meantime.mcmc
), 2)
df.estimator <- round(data.frame(
  raw.gvam = c(
    unlist(res.gvam$beta),
    unlist(res.gvam$sigma),
    res.gvam$rho),
  raw.mcmc = c(
    unlist(res.mcmc$beta),
    unlist(res.mcmc$sigma),
    res.mcmc$rho),
  mean.gvam = c(
    unlist(mean.beta.gvam),
    unlist(mean.sigma.gvam),
    mean.rho.gvam),
  mean.mcmc = c(
    unlist(mean.beta.mcmc),
    unlist(mean.sigma.mcmc),
    mean.rho.mcmc),
  sd.gvam = c(
    unlist(sd.beta.gvam),
    unlist(sd.sigma.gvam),
    sd.rho.gvam),
  sd.mcmc = c(
    unlist(sd.beta.mcmc),
    unlist(sd.sigma.mcmc),
    sd.rho.mcmc),
  conf2.5.gvam = c(
    unlist(conf.beta.gvam)[names(unlist(conf.beta.gvam))=='2.5%'],
    unlist(conf.sigma.gvam)[names(unlist(conf.sigma.gvam))=='2.5%'],
    conf.rho.gvam[names(conf.rho.gvam)=='2.5%']),
  conf97.5.gvam = c(
    unlist(conf.beta.gvam)[names(unlist(conf.beta.gvam))=='97.5%'],
    unlist(conf.sigma.gvam)[names(unlist(conf.sigma.gvam))=='97.5%'],
    conf.rho.gvam[names(conf.rho.gvam)=='97.5%']),
  conf2.5.mcmc = c(
    unlist(conf.beta.mcmc)[names(unlist(conf.beta.mcmc))=='2.5%'],
    unlist(conf.sigma.mcmc)[names(unlist(conf.sigma.mcmc))=='2.5%'],
    conf.rho.mcmc[names(conf.rho.mcmc)=='2.5%']),
  conf97.5.mcmc = c(
    unlist(conf.beta.mcmc)[names(unlist(conf.beta.mcmc))=='97.5%'],
    unlist(conf.sigma.mcmc)[names(unlist(conf.sigma.mcmc))=='97.5%'],
    conf.rho.mcmc[names(conf.rho.mcmc)=='97.5%'])
), 5)
rownames(df.estimator) <- c(
  paste0('beta1_',0:(ncol(fixed.MM[[1]])-1)),
  paste0('beta2_',0:(ncol(fixed.MM[[2]])-1)),
  'sigma1^2', 'sigma2^2', 'rho')

print(df.time)
print(df.estimator)

##########




