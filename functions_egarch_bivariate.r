#-------------------------------------------------------------------------------
#GARCH log-likelihood calculation
#-------------------------------------------------------------------------------

#NOTE: for NGARCH mean, specify mean_incl=FALSE and model$mean="r+u" as this implies u=0
EGARCH_bivariate <- function(all_vars,
                             mean_incl=TRUE,
                             model=list('vol'="EGARCH",
                                        'mean'=c("u", "r+u"),
                                        'lambda'=c("none", "NGARCH")
                                        ),
                             init=list(method=c('estim', 'unc', 'fixed'),
                                       h_0=h_0
                                       ),
                             logr=0,
                             div=0,
                             logR,
                             mult=100,
                             output=c("loglik","full")){

  if(model$mean=="r+u"){
    logR <- logR + div - logr #excess returns
  } 
  #note: when model$mean="u", we ignore dividends and interest rates

  E <- 2/sqrt(2*pi) #E[|z|]
    
  params <- working_to_real_EGARCH(all_vars, mean_incl, model, init)
  u  <- params$u; a0 <- params$a0; a1 <- params$a1; b <- params$b;
  g1 <- params$g1; rho <- params$rho; lambda <- params$lambda; 
  h_0 <- params$h_0; 
  
  h_t <- h_0 #h_0 represents conditional variance at t=1
  
  if(output=="full"){
    mu_all <- matrix(ncol=2, nrow=nrow(logR))
    colnames(mu_all) <- colnames(logR)
    rownames(mu_all) <- rownames(logR) 
    
    h_all <- matrix(ncol=2, nrow=nrow(logR)+1)
    colnames(h_all) <- colnames(logR)
    rownames(h_all) <- c(rownames(logR), "oos")
    h_all[1,] <- h_t
  }
   
  loglik <- 0
  for(i in 1:nrow(logR)) {
    #conditional mean
    if(model$lambda=='NGARCH') {
      mu_t <- u + lambda*sqrt(h_t) - 1/2*h_t/mult
      #division of sig_2 by mult to preserve consistency
    } else if(model$lambda=='none'){
      mu_t <- u
    }
    
    #bivariate log-density
    lognorm2 <- -log(2*pi) - 0.5*log(h_t[1]) - 0.5*log(h_t[2]) - 0.5*log(1-rho^2) -  
                1/(2*(1-rho^2)) * ( (logR[i,1]-mu_t[1])^2/h_t[1] + 
                                    (logR[i,2]-mu_t[2])^2/h_t[2] -
                                     2*rho*(logR[i,1]-mu_t[1])*(logR[i,2]-mu_t[2])/sqrt(h_t[1]*h_t[2]) )   
    
    #conditional covariance matrix
    #varcov <- diag(h_t)
    #varcov[1,2] <- varcov[2,1] <- rho*sqrt(h_t[1]*h_t[2])
    #lognorm2 <- dmvnorm(logR[i,], mean=mu_t, sigma=varcov, log=TRUE)
    
    #loglik
    loglik <- loglik + as.vector(lognorm2) 
    
    #vol update step    
    z <- (logR[i,]-mu_t) / sqrt(h_t)    
    h_t <- exp(a0 + a1*z + g1*(abs(z) - E) + b*log(h_t))
  
    if(output=="full"){
      mu_all[i,] <- mu_t
      h_all[i+1,] <- h_t
    }
  }

  if(is.na(loglik)) loglik <- -Inf 
  
  if(output=="loglik"){
    return(-loglik)
  } else if(output=="full"){
    return(list(loglik=-loglik, h_t=h_all, mu_t=mu_all))
  } 

}

#-------------------------------------------------------------------------------
#FUNCTIONS TO CONVERT PARAMETERS
#-------------------------------------------------------------------------------

real_to_working_EGARCH <- function(u, a0, a1, b, g1=0, lambda=0, rho=0,
                                   mean_incl=TRUE, model, init){

  if(mean_incl==FALSE) u <- numeric(0)
  
  if(model$lambda=="none") lambda <- numeric(0)
  
  b.tr <- sqrt(1/b - 1) #[0,1]

  rho.tr <- tan(rho*pi/2) #rho in [-1,1]

  if(init$method=='estim') h_0.tr <- sqrt(init$h_0) else h_0.tr <- numeric(0)

  all_vars <- c(u, a0, a1, b.tr, g1, lambda, rho.tr, h_0.tr)

  return(as.vector(all_vars))

}

working_to_real_EGARCH <- function(all_vars, mean_incl=TRUE, model, init){

  if(mean_incl==TRUE){
    u <- all_vars[1:2]
    all_vars <- all_vars[-c(1:2)]
  } else u <- rep(0,2)

  a0 <- all_vars[1:2]
  all_vars <- all_vars[-c(1:2)]
  
  a1 <- all_vars[1:2]
  all_vars <- all_vars[-c(1:2)]

  b <- 1 / (1+all_vars[1:2]^2)
  all_vars <- all_vars[-c(1:2)]
  
  g1 <- all_vars[1:2]
  all_vars <- all_vars[-c(1:2)]

  if(model$lambda=="NGARCH") {
    lambda <- all_vars[1:2]
    all_vars <- all_vars[-c(1:2)]
  } else if(model$lambda == "none") {
    lambda <- rep(0,2)  
  } 
 
  rho <- atan(all_vars[1])*2/pi
  all_vars <- all_vars[-1]
 
  #h_0 represents conditional variance at t=1
  if(init$method=='estim') {
    h_0 <- all_vars^2
  } else if(init$method=='unc'){
    h_0 <- exp( a0/(1-b) )
    #note: this is not equal to unconditional variance
  } else if(init$method=='fixed'){
    h_0 <- as.vector(init$h_0)
  }
    
  return(list(u=u, a0=a0, a1=a1, b=b, g1=g1, lambda=lambda, rho=rho, h_0=h_0))
}

#jacobian of the working_to_real transformation
working_to_real_EGARCH_jacobian <- function(u, a0, a1, b, g1=0, lambda=0, rho=0,
                                            mean_incl=TRUE, model, init){
  
  if(mean_incl==TRUE) {
    jac <- rep(1,6)  #u, a0, a1
  } else {
    jac <- rep(1,4)  #a0, a1
  }  
  
  jac_b <- -2*b^(3/2)*sqrt(1-b)  
  jac <- c(jac, jac_b)
    
  if(model$lambda=="NGARCH") {
    jac <- c(jac, rep(1,4)) #g1, lambda
  } else if(model$lambda == "none") {
    jac <- c(jac, rep(1,2)) #g1
  } 
  
  jac_rho <- 2/pi * 1 / (1+(tan(rho*pi/2))^2)  
  jac <- c(jac, jac_rho)
  
  if(init$method=='estim') jac <- c(jac,2*sqrt(h_0))
  
  return(jac)
}

