
library(statmod)
source("CalculateB.Rs")

###############################################################################

GVA.randint.f <- function(vtheta,family,vy,mX,vn,id,gh) 
{
    p <- ncol(mX)
    m <- length(vn)
    
    vbeta <- vtheta[1:p]
    gamma <- vtheta[1+p]
    vxi   <- vtheta[1+p+(1:(2*m))]
    vmu   <- vxi[2*(1:m) - 1]
    vzeta <- vxi[2*(1:m)]
    vlambda <- exp(vzeta)
    
    veta.vmu    <- mX%*%vbeta + rep(vmu,vn)     
    veta.sigma2 <- rep(vlambda,vn)
    vB0 <- B0.fun(family,veta.vmu,veta.sigma2,gh) 
     
    sl <- sum(vy*veta.vmu - vB0) - 0.5*m*gamma + 0.5*sum(vzeta - exp(-gamma)*(vmu^2 + vlambda)) + m/2
    if (family=="POISSON") {
        sl <- sl - sum(lgamma(vy+1))
    }
    
    return(-sl)
}

###############################################################################

GVA.randint.NRDIR <- function(vtheta,family,vy,mX,vn,id,gh,linds) 
{
    p <- ncol(mX); m <- length(vn)
    
    vbeta <- vtheta[1:p]
    gamma <- vtheta[1+p]
    vxi   <- vtheta[1+p+(1:(2*m))]
    vmu   <- vxi[2*(1:m) - 1]
    vzeta <- vxi[2*(1:m)]
    
    sigma2  <- exp(gamma)
    vlambda <- exp(vzeta)
    
    EPS <- 1.0E-15
    
    ###########################################################################
    
    # Calculate 1d integrals and some auxiliary information

    veta.vmu    <- mX%*%vbeta + rep(vmu,vn)     
    veta.sigma2 <- rep(vlambda,vn)
    res  <- B1234.fun(family,veta.vmu,veta.sigma2,gh) 
    
    veps <- vy - res$vB1
    vw2  <- res$vB2; dim(vw2) <- NULL    
    vw3  <- res$vB3; dim(vw3) <- NULL        
    XW2  <- mX*vw2
    XW3  <- mX*vw3    
    
    sEu <- sum(vmu^2 + vlambda)
    vs2 <- matrix(0,m,1)
    vs3 <- matrix(0,m,1)
    vs4 <- matrix(0,m,1)
    vseps <- matrix(0,m,1)
    for (i in 1:m) {
				ind <- linds[[i]]
		    vs2[i] <-  sum(res$vB2[ind])
		    vs3[i] <-  sum(res$vB3[ind])
		    vs4[i] <-  sum(res$vB4[ind])
		    vseps[i]  <-  sum(veps[ind])
    }
    		    
    ###########################################################################
    	
    # Calculate Gradient Information
		vg.vbeta <- t(mX)%*%veps
		vg.gamma <- -0.5*m + 0.5*sEu/sigma2
		vg.vmu   <- vseps - vmu/sigma2
		vg.vzeta <- 0.5*(1  - vlambda*(1/sigma2 + vs2))		
		
		vg.vtheta <- matrix(c(vg.vbeta,vg.gamma))
		vg.vxi <- matrix(0,2*m,1)
    vg.vxi[2*(1:m)-1] <- vg.vmu
    vg.vxi[2*(1:m)]   <- vg.vzeta		
    
    vg <- matrix(c(vg.vtheta,vg.vxi))
        
    ###########################################################################
          
    # Allocate some memory
		mH.00 <- matrix(0,p+1,p+1)
		mH.0i <- matrix(0,p+1,2*m)
		mH.ii.inv <- matrix(0,2,2)			

		dvtheta <- matrix(0,p+1+2*m,1)   		
		mA      <- matrix(0,2*m,p+1)	 	
    
    # Calculate Search Direction
    
		mH.00[1:p,1:p] <-  -t(XW2)%*%mX - diag(EPS,p)
		mH.00[1+p,1+p] <-  -0.5*sEu/sigma2 - EPS
		
		for (i in 1:p) {
		    mH.0i[i,2*(1:m)-1] <- -tapply(XW2[,i],id,sum)
		    mH.0i[i,2*(1:m)]   <- -0.5*vlambda*tapply(XW3[,i],id,sum)
		}
		mH.0i[1+p,2*(1:m)-1] <- vmu/sigma2
		mH.0i[1+p,2*(1:m)]   <- 0.5*vlambda/sigma2
		
		mH.ii.11 <- -(1/sigma2) - vs2 - rep(EPS,m)
		mH.ii.22 <- -0.5*vlambda/sigma2 - 0.5*vlambda*(vs2 + 0.5*vlambda*vs4) - rep(EPS,m)
		mH.ii.12 <- -0.5*vlambda*vs3
		
		det.term <- mH.ii.11*mH.ii.22 - mH.ii.12^2
		
		vg.til <- vg.vtheta
		mH.til <- mH.00
		for (i in 1:m) {
		    
		    mH.ii.inv[1,1] <-  mH.ii.22[i]/det.term[i]
		    mH.ii.inv[2,2] <-  mH.ii.11[i]/det.term[i]
		    mH.ii.inv[1,2] <- -mH.ii.12[i]/det.term[i]
		    mH.ii.inv[2,1] <- mH.ii.inv[1,2]
		    
		    vg.i <- vg.vxi[2*(i-1)+(1:2)]
		    mH.i <- mH.0i[,2*(i-1)+(1:2)]
		    
		    vb <- mH.ii.inv%*%vg.i
		    mB <- mH.ii.inv%*%t(mH.i)
		    
		    dvtheta[1+p+2*(i-1)+(1:2)] <- vb
		    mA[2*(i-1)+(1:2),1:(1+p)]  <- mB
		    
		    vg.til <- vg.til - mH.i%*%vb
		    mH.til <- mH.til - mH.i%*%mB
		}
		
		sI <- solve(mH.til)		#ginv
		vs <- sI%*%vg.til
		dvtheta[1:(p+1)] <- vs 
		dvtheta[p+1+(1:(2*m))] <- dvtheta[p+1+(1:(2*m))] - mA%*%vs
							
    return(list(dvtheta=-dvtheta,vg=vg,sI=-sI))
}

###############################################################################

GVA.randint.FIT <- function(id,vy,mX,vbeta,sigma2,vmu,vlambda,
                             family = "POISSON",
                             MAXITER = 100,
                             EPS.TERM.PAR  = 1.0E-5,
                             EPS.TERM.GRAD = 1.0E-5) 
{    
    
    if (family=="LOGISTIC") {
				# Generate Quadrature Points
				N <- 15
				gh.temp <- gauss.quad(N,kind="hermite")
				gh <- list(x=gh.temp$nodes,w=gh.temp$weights,w.til=gh.temp$weights*exp(gh.temp$nodes^2))  
		}  

    p <- ncol(mX)
    unique.ids <- unique(id)
    m  <- length(unique.ids)
    vn <- as.vector(table(id))

    gamma <- log(sigma2)
    vzeta <- log(vlambda)
		vxi <- matrix(0,2*m,1)
    vxi[2*(1:m)-1] <- vmu
    vxi[2*(1:m)]   <- vzeta   
    
    
    linds <- list()
    for (i in 1:m) {
        linds[[i]]  <- which(id==unique.ids[i])
    }      
    
    vtheta <- c(vbeta,gamma,vxi)
    sl <- -GVA.randint.f(vtheta,family,vy,mX,vn,id,gh) 
    
    cat("Initial variational log-likelihood = ",sl)
    cat("\n\n")
    
    cat("ITER    VARLOGLIK        STEPSIZE         GRADERR       MAXPARDIF \n")
       
    for (ITER in 1:MAXITER) {
        vtheta.old <- vtheta
    
        # Calculate Newton-Raphson Direction 
        res <- GVA.randint.NRDIR(vtheta,family,vy,mX,vn,id,gh,linds)    
        
				stepsize <- 1                
        vtheta.temp <- vtheta + stepsize*res$dvtheta
        sl.temp <- -GVA.randint.f(vtheta.temp,family,vy,mX,vn,id,gh) 
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
		            vsl[JTER] <- -GVA.randint.f(vtheta.temp,family,vy,mX,vn,id,gh) 
		        }                   
		        sl <- max(vsl)
		        stepsize <- vstepsize[which.max(vsl)]
		        vtheta <- vtheta + stepsize*res$dvtheta
		    } else {
		        sl <- sl.temp
		        vtheta <- vtheta.temp           
        }
    
    
        # Check termination condtions
        err.par <- max( abs(vtheta - vtheta.old) / pmax(abs(vtheta.old),1e-10) )
        err.vg  <- sqrt(sum(res$vg^2)/length(vtheta))
        
        cat(sprintf("%4d    %E    %+E    %E  %E\n",ITER,sl,stepsize,err.vg,err.par))

        if (err.vg<EPS.TERM.GRAD) { break; }
        if (err.par<EPS.TERM.PAR) { break; }
        
        if (exp(vtheta[1+p])<5.0E-4) { break; } # Break if sigma2 is getting too small

    }
    
    vbeta <- vtheta[1:p]
    gamma <- vtheta[1+p]
    vmu   <- vtheta[p+2*(1:m)]
    vzeta <- vtheta[1+p+2*(1:m)]
    
		vxi <- matrix(0,2*m,1)
    vxi[2*(1:m)-1] <- vmu
    vxi[2*(1:m)]   <- exp(vzeta)   
    
    return(list(vbeta=vbeta,sigma2=exp(gamma),vmu=vmu,vlambda=exp(vzeta),vbeta.serr=sqrt(diag(res$sI[1:p,1:p])),elbo = sl,
    iterations = ITER))
}

###############################################################################