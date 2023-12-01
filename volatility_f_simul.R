
params <- working_to_real_EGARCH(Optim$pars, mean_incl, model, init)

params$a0[1] = -1.03122
params$a1[1] = -.09390
params$b[1] = 0.89453
params$lambda[1] = 0.06653
params$rho = 0
params$h_0[1] = 0.007529668

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
vol_final <- as.vector(sqrt(h_final*12))


# volatility generation ---------------------------------------------------
h <- c()
vol <- c()
omega <- -1.03122; alpha <- -0.09390; gamma <- 0.06653; beta <- 0.89453; h[1] <- sigma^2

for(j in 1:10000){
for(t in 2:312){
h[t] <- exp(omega + alpha * rnorm(1) + gamma * (abs(rnorm(1))-2/(sqrt(2*pi))) + beta* log(h[t-1]))}
vol[j] <- sqrt(12*h[length(h)])
}


2.59



