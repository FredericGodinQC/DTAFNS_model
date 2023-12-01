#-------------------------------------------------------------------------------
#LOAD FUNCTIONS AND DATA
#-------------------------------------------------------------------------------

#Install the following package
library(Rsolnp)

#Set path 
#path <- 'C:/Dropbox/MaciejEmmanuel/Rcode/WEBSITE/2_Equity_indices_EGARCH'
#path <- 'D:/Concordia Courses/Research paper/Research/Second/R_codes_mixed_bond/2_Equity_indices_EGARCH'
path <- 'D:/Concordia Courses/Research paper/Research/Second/Codes_myself/2_Equity_indices_EGARCH_mymodel'
setwd(path)

#Load functions
source("functions_egarch_bivariate.r")
source("expr_egarch_bivariate_load_data.r")
eval(expr_load_data_egarch)

#-------------------------------------------------------------------------------
#ESTIMATION OF EGARCH BIVARIATE
#-------------------------------------------------------------------------------

#-----------
#model entry
#-----------
mean_incl <- FALSE
model <- list('vol'="EGARCH", 'mean'="r+u", 'lambda'="NGARCH")
init <- list(method='fixed', h_0=diag(var(logR)))
init_par <- 1 #initial parameter set, see below

rel_tol <- 1e-8 #optimization parameter for solnp solver
mult <- 1       #logreturns are multiplied by mult before estimation

#--------------------
#load risk-free rates
#--------------------
logr <- ShortRateG3
logr <- logr[-length(logr)] / 12 #monthly rates
div <- 0 #dividend rate, assumed zero here

#----------------------
#load intial parameters
#----------------------if (mean_incl) u <- colMeans(logR) else u <- rep(0,2)
h_0 <- diag(var(logR))
if (init_par == 1) {
  a0 <- rep(-1, 2)
  a1 <- rep(-0.1, 2)
  b  <- rep(0.90, 2)
  g1 <- rep(0.30, 2)
  lambda <- rep(0.10, 2)
  rho <- 0.80
}

mean_ <- model$mean
lambda_ <- model$lambda
init_method <- init$method
model <- list('vol' = 'EGARCH', 'mean' = mean_, 'lambda' = lambda_)
init <- list(method = init_method, h_0 = h_0)

all_vars <- real_to_working_EGARCH(u, a0, a1, b, g1, lambda, rho, mean_incl, model, init)
loglik0 <- -EGARCH_bivariate(all_vars, mean_incl, model, init, 
                             logr=logr, div=div, logR=logR, mult=mult, output="loglik")

#----------------------
#estimation
#----------------------
start.time <- proc.time()["elapsed"]
Optim <- solnp(pars=all_vars, fun=EGARCH_bivariate, mean_incl=mean_incl,
                model=model, init=init, logr=logr, logR=logR, 
                  div=div, mult=mult, output="loglik", control=list(trace=0, tol=rel_tol))
end.time <- proc.time()["elapsed"]
run.time <- as.vector(end.time - start.time)

#----------------------
#results
#----------------------
params <- working_to_real_EGARCH(Optim$pars, mean_incl, model, init)
params_n <- length(Optim$pars)
u  <- params$u; a0 <- params$a0; a1 <- params$a1; b <- params$b;
g1 <- params$g1; rho <- params$rho; lambda <- params$lambda; 
h_0 <- params$h_0; 

all_vars <- real_to_working_EGARCH(u, a0, a1, b, g1, lambda, rho, mean_incl, model, 
                                    init=list(method = "fixed", h_0 = h_0))
output <- EGARCH_bivariate(all_vars, mean_incl, model, init=list(method = "fixed", h_0 = h_0), 
                           logr=logr, div=div, logR=logR, mult=mult, output="full")
loglik <- -output$loglik

h_final <- as.vector(output$h_t[rownames(logR)[nrow(logR)],]) #final monthly EGARCH variance
vol_final <- as.vector(sqrt(h_final*12)) #final annualized EGARCH volatility

res <- c('mean'=mean_,
         'u_estim'=mean_incl,
         'risk_premium'=lambda_,
         #'mult'=mult,
         'vol_init_method'=init$method,
         'init_par'=init_par,
         'loglik'=round(loglik, 3),         
         'AIC'=round(2*params_n - 2 * loglik, 3),
         'BIC'=round(log(nrow(logR))*params_n - 2 * loglik, 3),
         'par'=params_n,
          unlist(list(u=params$u, a0=params$a0, a1=params$a1, b=params$b,
                      g1=params$g1, lambda=params$lambda, rho=params$rho)),
         'vol_init'=round(sqrt(params$h_0*12),4),
         'vol_final'=round(vol_final,4),
         'h_final'=h_final,
         'n'=nrow(logR),
         #'loglik1'=round(loglik1,3),
         #'loglik100'=round(loglik100,3),
         'time'=round(run.time,2),
         'conv'=Optim$convergence
         )
as.data.frame(res)

params$h_0 <- h_final # I add this part to replaced h_0 with last h_t
save(params, file = paste0(getwd(),"/stock_param.RData"))

filename <- paste(colnames(logR), collapse="_")
filename <- paste(filename, "_", rownames(logR)[1], "_", 
                    rownames(logR)[nrow(logR)], "_", "EGARCH_single", sep="")
write.csv(res, paste(filename,"csv",sep="."))

#------------------------------------------------------
#standard errors
#------------------------------------------------------

#observed information matrix (working parameters)
I_obs <- Optim$hessian
#note: function returns negative of log-likelihood so
#we do not need to multiply by -1 to get information matrix

#asymptotic variance-covariance matrix (working parameters)
varcov_working <- solve(I_obs)
#jacobian matrix of transformation from working to real parameters
jac <- working_to_real_EGARCH_jacobian(u, a0, a1, b, g1, lambda, rho,
                                                   mean_incl, model, init)
#asymptotic variance-covariance matrix (real parameters)
varcov_real <- diag(jac) %*% varcov_working %*% diag(jac)
#standard errors
standard_errors <- sqrt(diag(varcov_real))
standard_errors_pars <- c("omega1", "omega2", "alpha1", "alpha2",
                          "beta1", "beta2", "gamma1", "gamma2",
                          "lambda1", "lambda2", "rho")
names(standard_errors) <- standard_errors_pars
as.matrix(standard_errors,ncol=1)
# omega1  0.55609040      
# omega2  0.56204440
# alpha1  0.06172138
# alpha2  0.07692681
# beta1   0.08539292
# beta2   0.08688848
# gamma1  0.07536069
# gamma2  0.07388792
# lambda1 0.03972168
# lambda2 0.04326677
# rho     0.02079423

#------------------------------------------------------
#outputs: log-returns, risk-free rates and volatilities
#------------------------------------------------------

output_vars <- cbind('t'=0:nrow(logR), 
                     'Rt_TSX'=c(NA,logR[,"TSX"]),  
                     'ht_TSX'=output$h_t[,"TSX"],
                     'volt_TSX'=sqrt(output$h_t[,"TSX"]*12),
                     'Rt_SP500'=c(NA,logR[,"SP500"]), 
                     'ht_SP500'=output$h_t[,"SP500"],
                     'volt_SP500'=sqrt(output$h_t[,"SP500"]*12),
                     'riskfree_annual'=ShortRateG3, 
                     'riskfree_monthly'=ShortRateG3/12)
rownames(output_vars) <- names(ShortRateG3)
filename <- paste(filename, "_", "output_vars", sep="")
write.csv(output_vars, paste(filename,"csv",sep="."))



