rm(list = ls())
getwd()
setwd("C:/Users/admin/Desktop")
library(psych)
library(dplyr)
library(Matrix)
library(MCMCglmm)
library(MASS)
library(haven)
library(glmmTMB)
source("CalculateB.Rs")
source("GVArandint.Rs")
source("mGLMM.Rs")
########### function ##########
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
########### get data ################
df_bac <- read.csv("taxonomy_counts_g.csv")
df_nut <- read.table("nutrition_totals.txt", header = TRUE, sep = "")

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

########### pre-disposal ##################
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

############# select #################
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


id <- as.vector(data$UserName)
Y <- as.matrix(data[c("y1","y2")])
fixed.MM <- list(
  as.matrix(cbind(1, data[c("S040","S080","S120","VITE_ADD")])),
  as.matrix(cbind(1, data[c('NIAC','P204')]))
)

############## fitting ##################
Xdependent <- FALSE  
beta0 <- 0  
gamma0 <- 0  
mu0   <- 0  
zeta0 <- 0  
family = "POISSON"

res.gvam <- GVAM.FIT(Xdependent=Xdependent, id=id, Y=Y, fixed.MM=fixed.MM, beta0=beta0, gamma0=gamma0, mu0=mu0, zeta0=zeta0, family = family)


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

########### result ############
describe(df_merge$y1) 
describe(df_merge$y2) 
describe(df_merge[c("S040","S080","S120","VITE_ADD",'NIAC','P204')])
print(res.gvam)
print(res.mcmc)
print(summary(model.mcmc))

