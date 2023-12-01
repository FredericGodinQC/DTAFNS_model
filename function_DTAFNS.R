
package_list <- c("dmm", "lubridate", "FKF", "matrixcalc", "Rsolnp", "DEoptim", "xtable")
for (i in package_list){
  if (require(i, character.only = TRUE) == FALSE)
  {
    install.packages(i);
    require(i, character.only = TRUE)
  }
}
library(dmm)
library(lubridate)
library(FKF)
library(matrixcalc)
library(Rsolnp)
library(DEoptim)
library(xtable)
lang <- "ENG"
constraints <- list(lambda1=0, 
                    lambda2="include", 
                    rho="include",
                    p=3)

unfactor_Data_Frame <- function(data_Frame, nb_Col_Begin = 1)
{
  library(dmm)
  nb_Cols <- ncol(data_Frame)
  for(i in nb_Col_Begin : nb_Cols)
  {
    if(class(data_Frame[, i]) == "factor")
    {
      data_Frame[, i] <- dmm::unfactor(data_Frame[, i])
    }
  }
  
  return(data_Frame)
}

extraire_Information_Annees_Mois_Jours_Date <- function(vecteur_Dates)
{
  library(lubridate)
  vecteur_Dates <- as.Date(vecteur_Dates)
  weekdays <- wday(vecteur_Dates)
  yeardays <- yday(vecteur_Dates)
  monthdays <- mday(vecteur_Dates)
  mois <- months(vecteur_Dates)
  trimestre <- quarter(vecteur_Dates)
  annees <- year(vecteur_Dates)
  return(list(weekdays = weekdays, yeardays = yeardays, monthdays = monthdays,
              mois = mois, trimestre = trimestre, annees = annees))
}

objet_Info_Date_Plus_Const <- function(objet_Info_Date, const = 15, lang="ENG")
{
  annees <- objet_Info_Date$annees
  if(lang=="ENG"){
    mois_lang <- c("January", "February", "March", "April", "May", "June",
                   "July", "August", "September", "October", "November",
                   "December")
  } else if(lang=="FRA"){
    mois_lang <- c("janvier", paste0("f", '\U00E9', "vrier"), "mars",
                   "avril", "mai", "juin",
                   "juillet", paste0("ao", '\U00FB', "t"), "septembre",
                   "octobre", "novembre", paste0("d", '\U00E9', "cembre"))
  }
  
  mois <- factor(objet_Info_Date$mois,
                 levels = mois_lang,
                 labels = 1 : 12)
  jours_Mois <- objet_Info_Date$monthdays
  
  matrice_Premier_Jours <- plus_Proche_15_Date_Mois(mois = mois, jours_Mois = jours_Mois, annees = annees, const = const)
  
  date_Premier_Mois <- as.Date(paste(matrice_Premier_Jours$annees, matrice_Premier_Jours$mois,
                                     matrice_Premier_Jours$jours_Mois, sep = "-"), format = "%Y-%m-%d")
  
  return(date_Premier_Mois)
}

plus_Proche_15_Date_Mois <- function(mois, jours_Mois, annees, const = 15)
{
  matrice_Donnees <- data.frame(annees, mois, jours_Mois)
  aggregate(data = matrice_Donnees, jours_Mois ~ mois + annees, FUN = plus_Proche_Const, const = const)
}

plus_Proche_Const <- function(vecteur, const = 15)
{
  vecteur[which.min(abs(vecteur - const))]
}

load_data2 <- expression({
  #options(digits = 13)
  yield_curve_Can <- read.table("yield_curves.csv", sep = ",", 
                                header = TRUE)
  dates <- as.Date(yield_curve_Can$Date)
  taux_Canada <- unfactor_Data_Frame(yield_curve_Can[, -c(1, 
                                                          ncol(yield_curve_Can))])
  ind <- which(dates == start.date):which(dates == end.date)
  dates <- dates[ind]
  taux_Canada <- taux_Canada[ind, ]
  taux_Canada[taux_Canada == " na"] <- NA 
  row_remove <- which(rowSums(is.na(taux_Canada)) == ncol(taux_Canada))
  dates <- dates[-row_remove]
  taux_Canada <- taux_Canada[-row_remove, ]
  #taux_Canada <- na.omit(taux_Canada)
  objet_Info_Date <- extraire_Information_Annees_Mois_Jours_Date(dates)
  date_Milieu_Mois <- objet_Info_Date_Plus_Const(objet_Info_Date, 
                                                 const = 31, lang = lang)
  #x <- seq(1,429,by = 12)
  #date_Milieu_Mois <- date_Milieu_Mois[x]   #Selecting yearly data
  index_Bool_Date <- as.Date(dates) %in% date_Milieu_Mois
  donnees_Calibration <- taux_Canada[index_Bool_Date, ]
  maturity_all <- seq(from = 0.25, to = 30, by = 0.25)
  maturity_all_years <- c(0.25, 0.5, 0.75, seq(from = 1, to = 30,
                                               by = 1)) 
  
  maturity_BoC <- c(0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 
                    30)
  if (maturity_type == "all") {
    maturity <- maturity_all
  }else if (maturity_type == "all_years") {
    maturity <- maturity_all_years
  }else if (maturity_type == "BoC") {
    maturity <- maturity_BoC
  }
  index <- which(maturity_all %in% maturity)
  y <- t(donnees_Calibration[, index])
  y <- matrix(as.numeric(y), nr = length(maturity))
    maturity <- maturity/delta_Temps
})

zeta <- function(r,maturity,number){
  if(number == 0){
    zeta <- (r-r^maturity)/(1-r)
  }else if(number == 1){
    zeta <- (r-maturity*(r^maturity)+(maturity-1)*r^(maturity+1))/(1-r)^2
  }else{
    zeta <- (-(maturity-1)^2*r^(maturity+2)+(2*maturity^2-2*maturity-1)*r^(maturity+1)-maturity^2*r^maturity+r^2+r)/(1-r)^3
  }
  return(zeta)
}

upsilon_11 <- function(maturity, sigma_11){
  sigma_11^2*maturity*(maturity-1)*(2*maturity-1)/6
}
upsilon_22 <- function(lambda, maturity, sigma_22){
  (sigma_22^2/lambda^2)*(maturity-2*
                           ((1-(1-lambda)^maturity)/lambda)+(1-(1-lambda)^(2*maturity))/(1-(1-lambda)^2))
}
upsilon_33 <- function(lambda, maturity, sigma_33){
  (sigma_33^2/lambda^2)*(maturity-2+zeta(((1-lambda)^2),(maturity -1),0)+lambda^2*zeta((1-lambda)^2,(maturity - 1),2)-
                           2*zeta((1-lambda),(maturity - 1),0)-2*lambda*zeta((1-lambda),(maturity-1),1)+2*lambda*zeta((1-lambda)^2,(maturity-1),1))
}
upsilon_12 <- function(rho_12,sigma_11,sigma_22,lambda,maturity){
  rho_12*sigma_11*sigma_22/lambda*(maturity*(maturity-1)/2-zeta((1-lambda),maturity,1))
}
upsilon_13 <- function(rho_13,sigma_11,sigma_33,lambda,maturity){
  rho_13*sigma_11*sigma_33/lambda*(maturity*(maturity-1)/2-1- zeta((1-lambda),(maturity-1),0)
                                   - (lambda+1)*zeta((1-lambda),(maturity-1),1) -
                                     lambda*zeta((1-lambda),(maturity-1),2))
}
upsilon_23 <- function(rho_23,sigma_22,sigma_33,lambda,maturity){
  rho_23*sigma_22*sigma_33/lambda*((maturity-2-(2-lambda)*zeta((1-lambda),(maturity-1),0) + 
                                      (1-lambda)*zeta((1-lambda)^2,(maturity-1),0))/lambda -
                                     zeta((1-lambda),(maturity-1),1) + (1-lambda)*zeta((1-lambda)^2,(maturity-1),1))
}  

parameter_load <- expression({
  rho <- param$rho
  rho_fin <- matrix(c(1,rho[1],rho[3],rho[1],1,rho[2],rho[3],rho[2],1),nr = 3, byrow = T)
  rho_12 <- rho_fin[1,2]; rho_13 <- rho_fin[1,3]; rho_23 <- rho_fin[2,3]
  theta_q <- param$theta_q
  sigma <- diag(param$sigma)
  sigma_11 <- sigma[1,1]; sigma_22 <- sigma[2,2]; sigma_33 <- sigma[3,3]
  lambda <- param$lambda
  
  eta <- diag(param$eta)
  h <- param$h
  if(lambda == (0)){lambda = 1}
  kappa_q <- matrix(c(0,0,0,0,lambda,-lambda,0,0,lambda), nr = 3, byrow = T)
  kappa_p <- kappa_q + sigma %*% eta
})


construct_model_parameter <- function(param, maturity) 
{
  p = 3
  delta <- 1/12
  nb_Maturity <- length(maturity)
  eval(parameter_load)
  #state_mean <- theta_p
  #state_var <- sigma^2 %*% solve(1 - (1 - kappa_p)^2)
  #a0 <- state_mean
  #P0 <- state_var
  a0 <- param$a0
  P0 <- diag(param$P0)
  #lambda <- lambda* delta_Temps # add this line for test
  logAA <- -delta*theta_q[2]*(maturity - 1 - (1 - lambda - (1-lambda)^maturity)/lambda) - delta*theta_q[3]*
    ((maturity-1)*(1-lambda)^(maturity-1) - (1-(1- lambda)^(maturity - 1))/lambda) + 
    0.5*delta^2*(upsilon_11(maturity = maturity, sigma_11 = sigma_11) + upsilon_22(lambda = lambda, maturity = maturity, sigma_22 = sigma_22) + upsilon_33(lambda, maturity, sigma_33) +
                   2*(upsilon_12(rho_12 = rho_12,sigma_11 = sigma_11,sigma_22 = sigma_22,lambda = lambda,maturity = maturity) +
                        upsilon_13(rho_13 = rho_13,sigma_11 = sigma_11,sigma_33 = sigma_33,lambda = lambda,maturity = maturity) +
                        upsilon_23(rho_23 = rho_23,sigma_22 = sigma_22,sigma_33 = sigma_33,lambda = lambda,maturity = maturity)))
  #logAA <- 0 #test DNS
  
  BB_1 <- delta * maturity
  BB_2 <- delta*(1-(1-lambda)^maturity)/lambda
  BB_3 <- delta*((1-(1-lambda)^(maturity-1))/lambda - (maturity-1) * (1-lambda)^(maturity-1))
  #BB_2 <- delta*(1-exp(-lambda*maturity))/lambda    # test  DNS
  #BB_3 <- delta*(1-exp(-lambda*maturity))/lambda  - delta*maturity*exp(-lambda*maturity)#test DNS
  aa <- -logAA/(delta*maturity)
  BB <- matrix(c(BB_1/(delta*maturity), BB_2/(delta*maturity), BB_3/(delta*maturity)),nc = 3, byrow = F)
  ct <- aa
  Zt <- BB
  #h.tr <- sqrt(1/h - 16)
  GGt <- diag(rep(h, nb_Maturity))
  HHt <- t(sigma) %*% rho_fin %*% sigma
  Tt <- diag(3) - kappa_p
  #dt <- kappa_p %*% theta_p
  dt <- kappa_q %*% theta_q
  
  
  return(list(ct = ct, Zt = Zt, GGt = GGt, dt = dt, Tt = Tt, 
              HHt = HHt, a0 = a0, P0 = P0, maturity = maturity, delta_Temps = delta))
}


g3d_log_like_2 <- function (par, y, maturity, delta_Temps, constraints, 
                            a0 = NULL, P0 = NULL,par_input="working", output = "loglik", Type = "discrete") 
{
  
  if(par_input == "working"){
    param <- list(theta_q = par[1:3], sigma=par[4:6],
                  lambda = par[7],
                  rho = par[8:10],
                  eta =par[11:13],
                  h = par[14],
                  a0 = par[15:17],
                  P0 = par[18:20],
                  delta_Temps=par[21])
  }else{param <- par}
  p <- 3
  n <- ncol(y)
  d <- nrow(y)
  liste_Param_Espace_Etat <- construct_model_parameter(param = param,maturity = maturity) 
  
  dt <- liste_Param_Espace_Etat$dt
  Zt <- liste_Param_Espace_Etat$Zt
  ct <- liste_Param_Espace_Etat$ct
  Tt <- liste_Param_Espace_Etat$Tt
  HHt <- liste_Param_Espace_Etat$HHt
  GGt <- liste_Param_Espace_Etat$GGt
  if (is.null(a0)){ 
    a0 <- liste_Param_Espace_Etat$a0}else{a0 <- a0[1:p]}
  
  if(is.null(P0)){ 
    P0 <- liste_Param_Espace_Etat$P0}else{P0 <- as.matrix(as.matrix(P0)[1:p, 1:p])}
  
  
  m <- constraints$p
  eigen_Values <- try(eigen(HHt, only.values = TRUE), silent = TRUE)
  if (is.null(attr(eigen_Values, "class")) == FALSE) {
    return(10^10)
  }else {
    bool_Eigen_Neg <- any(eigen_Values$values < 0)
    if (bool_Eigen_Neg == TRUE) 
      return(10^10)
  }
  
  if ((length(dt) != m)) {
    stop("the dimension of dt is wrong")
  }else if ((dim(as.matrix(Tt))[1] != m) | (dim(as.matrix(Tt))[2] != 
                                            m)) {
    stop("the dimension of Tt is wrong")
  }else if ((dim(as.matrix(HHt))[1] != m) | (dim(as.matrix(HHt))[2] != 
                                             m)) {
    stop("the dimension of HHt is wrong")
  }else if ((length(ct) != d)) {
    stop("the dimension of ct is wrong")
  }else if ((dim(as.matrix(Zt))[1] != d) | (dim(as.matrix(Zt))[2] != 
                                            m)) {
    stop("the dimension of Zt is wrong")
  }else if ((dim(as.matrix(Zt))[1] != d) | (dim(as.matrix(Zt))[2] != 
                                            m)) {
    stop("the dimension of Zt is wrong")
  }else if ((dim(as.matrix(GGt))[1] != d) | (dim(as.matrix(GGt))[2] != 
                                             d)) {
    stop("the dimension of GGt is wrong")
  }else if ((dim(as.matrix(P0))[2] != m) | (dim(as.matrix(P0))[1] != 
                                            m)) {
    stop("the dimension of P0 is wrong")
  }
  
  out <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, 
             Zt = Zt, HHt = HHt, GGt = GGt, yt = y)
  if(output == "loglik"){
    logLik <- out$logLik
    if(is.na(logLik) == TRUE){
      return(10^10)
    }else{return(-logLik)} 
  }else{
    states_filter_mean <- out$att
    states_filter_var <- out$Ptt
    states_pred_mean <- out$at
    states_pred_var <- out$Pt
    error_pred_mean <- out$vt
    states_smooth_mean <- states_filter_mean
    states_smooth_var <- states_filter_var
    for (i in ncol(y):2) {
      J <- states_filter_var[, , i - 1] %*% t(Tt) %*% 
        solve(states_pred_var[, , i])
      states_smooth_mean[, i - 1] <- states_filter_mean[, 
                                                        i - 1] + J %*% (states_smooth_mean[, i] - states_pred_mean[, 
                                                                                                                   i])
      states_smooth_var[, , i - 1] <- states_filter_var[, 
                                                        , i - 1] + J %*% (states_smooth_var[, , i] - 
                                                                            states_pred_var[, , i]) %*% t(J)
    }
    
    states_filter_mean <- data.frame(t(states_filter_mean))
    dimnames(states_filter_mean)[[2]] <- paste0("x", 1:p, 
                                                "_filter")
    
   states_filter_var <- data.frame(t(apply(states_filter_var, 
                                              3, diag)))
    
    
    dimnames(states_filter_var)[[2]] <- paste0("x", 1:p, 
                                               "_filter_var")
    states_pred_mean <- data.frame(t(states_pred_mean))
    dimnames(states_pred_mean)[[2]] <- paste0("x", 1:p, 
                                              "_pred")
    
   states_pred_var <- data.frame(t(apply(states_pred_var, 
                                            3, diag)))
    
    
    dimnames(states_pred_var)[[2]] <- paste0("x", 1:p, "_pred_var")
    states_smooth_mean <- data.frame(t(states_smooth_mean))
    dimnames(states_smooth_mean)[[2]] <- paste0("x", 1:p, 
                                                "_smooth")
  
      states_smooth_var <- data.frame(t(apply(states_smooth_var, 
                                              3, diag)))
  
    
    dimnames(states_smooth_var)[[2]] <- paste0("x", 1:p, 
                                               "_smooth_var")
    ticker <- c(date_Milieu_Mois, NA)
    states_filter_mean <- rbind(states_filter_mean, NA)
    states_filter_var <- rbind(states_filter_var, NA)
    states_smooth_mean <- rbind(states_smooth_mean, NA)
    states_smooth_var <- rbind(states_smooth_var, NA)
    results <- data.frame(ticker, states_filter_mean, states_filter_var, 
                          states_pred_mean, states_pred_var, states_smooth_mean, 
                          states_smooth_var)
    out$summary <- results
    error_pred_mean <- data.frame(ticker[-length(ticker)], 
                                  t(error_pred_mean))
    dimnames(error_pred_mean)[[2]] <- c("date", paste0("maturity_", 
                                                       maturity_all_years))
    out$error_pred <- error_pred_mean
    return(out)
  }
}


param_work <- function(param){
  par <- unlist(param, use.names = FALSE)
  return(par)
}

work_param <- function(par){
  param <- list(theta_q = par[1:3], sigma=par[4:6],
                lambda = par[7],
                rho = par[8:10],
                eta = par[11:13], h = par[14],
                a0 = par[15:17],
                P0 = par[18:20],
                delta_Temps=par[21])
  return(param)
}



ineqn1 <- function(par, y, maturity, model_type, 
                   delta_Temps, constraints, a0=a0, P0=P0,
                   par_input="working",output="loglik") 
{
  param <- list(theta_q = par[1:3], sigma=par[4:6],
                lambda = par[7],
                rho = par[8:10],
                eta = par[11:13], h = par[14],
                a0 = par[15:17],
                P0 = par[18:20],
                delta_Temps=par[21])
  delta_Temps <- param$delta_Temps
  sigma <- param$sigma
  lambda <- param$lambda
  eta <- diag(param$eta)
  sigma <- diag(param$sigma)
  #sigma_11 <- sigma[1,1]; sigma_22 <- sigma[2,2]; sigma_33 <- sigma[3,3]
  rho <- param$rho
  #kappa_q <- matrix(c(0,0,0,0,lambda,-lambda,0,0,lambda), nr = 3, byrow = T)
  #kappa_p <- kappa_q + sigma %*% eta
  rho_fin <- matrix(c(1,rho[1],rho[3],rho[1],1,rho[2],rho[3],rho[2],1),nr = 3, byrow = T)
  #new version
  up11 <- upsilon_11(maturity, sigma[1,1])
  up22 <- upsilon_22(lambda, maturity, sigma_22 = sigma[2,2])
  up33 <- upsilon_33(lambda, maturity, sigma_33 = sigma[3,3])
  up12 <- upsilon_12(rho_fin[1,2],sigma_11 = sigma[1,1],sigma_22 = sigma[2,2],lambda,maturity)
  up13 <- upsilon_13(rho_fin[1,3],sigma_11 = sigma[1,1],sigma_33 = sigma[3,3],lambda,maturity)
  up23 <- upsilon_23(rho_fin[2,3],sigma_22 = sigma[2,2],sigma_33 = sigma[3,3],lambda,maturity)
  rou12 <- up12/(sqrt(abs(up11))*sqrt(abs(up22)))
  rou13 <- up13/(sqrt(abs(up11))*sqrt(abs(up33)))
  rou23 <-  up23/(sqrt(abs(up33))*sqrt(abs(up22)))
  
  #return(c(up33))
  return(c(rou12,rou13,rou23))
  
  # old version
  #sigma2_3 <- sigma2_3_f(sigma, maturity, lambda)
  #psi_23 <- psi_23_f(lambda, maturity, rho_fin, sigma)
  #sigma2_2 <- psi_12_f(lambda, maturity, sigma,rho_fin)
  #psi_12 <- psi_12_f(lambda, maturity, sigma,rho_fin)
  #return(c(sigma2_3,psi_23, sigma2_2, psi_12))
}



#ineqn1(par, y, maturity, model_type,  delta_Temps, constraints, a0=a0, P0=P0,  par_input="working",output="loglik") 



convert_q_to_p <- function(param){
  rho <- param$rho
  rho_fin <- matrix(c(1,rho[1],rho[3],rho[1],1,rho[2],rho[3],rho[2],1),nr = 3, byrow = T)
  theta_q <- param$theta_q
  sigma <- diag(param$sigma)
  eta <- diag(param$eta)
  h <- param$h
  lambda <- param$lambda
  kappa_q <- matrix(c(0,0,0,0,lambda,-lambda,0,0,lambda), nr = 3, byrow = T)
  kappa_p <- kappa_q + sigma %*% eta
  #theta_p <- solve(kappa_p)%*%kappa_q%*%theta_q
  theta_p <- c(0, (lambda*theta_q[2]-lambda*theta_q[3]*sigma[3,3]*eta[3,3]/kappa_p[3,3])/kappa_p[2,2],
               lambda*theta_q[3]/kappa_p[3,3])
  param_p <- list(rho = rho_fin,  kappa_p =  kappa_p, sigma = sigma,theta_p = theta_p, lambda = lambda, h = h, eta = eta)
  return(param_p)
}
out_smaple_loglike <- function(param = param, type, start = "2018-11-1"){
  
  start.date <- as.Date(start)  # "1986-01-01"
  end.date <- as.Date("2022-02-09")
  eval(load_data2)
  P0 <- diag(c(0.002^2, 0.002^2, 0.002^2))
  model_type <- "discrete"
  constraints <- list(lambda1=0, 
                      lambda2="include", 
                      rho="include",
                      p=3)
  if(type == "benchmark"){
    a0 <- c(0.04, 0.03, 0.03)  
    param_working <- real_to_working(param, constraints)
    loglik0 <- -g3d_Log_Lik(param_working, y, maturity, model_type, 
                            delta_Temps, constraints, a0=a0, P0=P0,
                            par_input="working", output="loglik")
  }else if(type == "Model"){
    a0 <- c(0.05, 0.04, 0.03)  
    par <- param_work(param)
    loglik0 <- -g3d_log_like_2(par, y, maturity, model_type, 
                               delta_Temps, constraints, a0=a0, P0=P0,
                               par_input="working",output="loglik")
  }
  return(loglik0)
}

bic <- function(logl, n, k){
  k*log(n) - 2* logl
}

# Out of sample mean error
out_sample_error <- function(param, type, start = "2018-11-1"){
  start.date <- as.Date(start)  # "1986-01-01"
  end.date <- as.Date("2022-02-09")
  eval(load_data2)
  P0 <- diag(c(0.002^2, 0.002^2, 0.002^2))
  model_type <- "discrete"
  constraints <- list(lambda1=0, 
                      lambda2="include", 
                      rho="include",
                      p=3)
  if(type == "benchmark"){
    # benchmark is godin paper model
    p=3
    a0 <- c(-0.06896636, -0.00619378, 0.09397286)  
    param <- working_to_real(Optim_bench$pars, model_type, delta_Temps, constraints)
    results <- g3d_Log_Lik(param, y, maturity, model_type, delta_Temps, constraints,
                           a0=a0, P0=P0,
                           par_input="real", output="full")
    output <- results$error_pred
  }else if(type == "Model"){
    a0 <- c(-0.02851708, 0.04856946, 0.06650171)  
    #par <- param_work(param)
    results <- g3d_log_like_2(param, y, maturity, model_type, 
                              delta_Temps, constraints, a0=a0,
                              par_input="real",output="full")
    output <- results$error_pred
  }else{
    a0 <- c(-3.704684, 3.724725, 3.977532)  
    #par <- param_work(param)
    results <- g3d_log_like_dns(par=dns_param, y, maturity, model_type, delta_Temps, constraints, 
                                a0 = a0, P0 = P0, output = "full", Type = "discrete") 
    output <- results$error_pred
  }
  
  mae <- colMeans(abs(output[,-1]))
  rmse <- sqrt(colMeans(output[,-1]^2))
  mean <- colMeans(output[,-1])
  
  
  return(list(Mean = mean,MAE = mae,RMSE = rmse))
}
table_error <- function(benchmark = benchmark_error, model = model_error){
  library(xtable)
  maturity <- c(3,6,9,seq(12,360,by=12))
  #table <- 100*cbind(Mean = as.numeric(model_error$Mean), MAE = as.numeric(model_error$MAE),
  #              RMSE = as.numeric(model_error$RMSE), Mean = as.numeric(benchmark_error$Mean), MAE = as.numeric(benchmark_error$MAE),
  #             RMSE = as.numeric(benchmark_error$RMSE),
  #            rMAE = as.numeric(model_error$MAE)/as.numeric(benchmark_error$MAE))
  table <- 100*cbind(RMSE = as.numeric(model$RMSE), MAE = as.numeric(model$MAE),
                     RMSE = as.numeric(benchmark$RMSE), MAE = as.numeric(benchmark$MAE),
                     rMAE = as.numeric(model$MAE)/as.numeric(benchmark$MAE))
  rownames(table) <- maturity
  maturity2 <- c(maturity[18:34])
  table2 <- cbind(table[1:17,],maturity2, table[c(18:33,NA), ])
  xtable(table2)
  
}

optim_nelson_sigels <- function(par, maturity, y){
  rates <- t(y)
  alpha_1 <- par[1]
  alpha_2 <- par[2]
  alpha_3 <- par[3]
  lambda <- par[4]
  implied_curve_ns <- alpha_1 + alpha_2*((1-exp(-lambda*maturity))/lambda*maturity) + 
    alpha_3*(exp(-lambda*maturity)/lambda*maturity-exp(-lambda*maturity)) 
  error <- sum(abs(rates - implied_curve_ns)/rates, na.rm = T)
  return(error)
}

nelson_sigels <- function(par, maturity){
  alpha_1 <- par[1]
  alpha_2 <- par[2]
  alpha_3 <- par[3]
  lambda <- par[4]
  implied_curve_ns <- alpha_1 + alpha_2*((1-exp(-lambda*maturity))/lambda*maturity) + 
    alpha_3*(exp(-lambda*maturity)/lambda*maturity-exp(-lambda*maturity)) 
  return(implied_curve_ns)
}

nelson_sigels_error <- function(par, model){
  start.date <- as.Date("2018-10-31")  # "1986-01-01"
  end.date <- as.Date("2022-02-09")
  eval(load_data2)
  rates <- t(y)
  model <- matrix(c(rep(model,NROW(rates))),byrow = T, nr = NROW(rates))
  error <- rates - model
  mae <- colMeans(abs(error),na.rm = T)
  rmse <- sqrt(colMeans(error^2, na.rm = T))
  mean <- colMeans(error,na.rm = T)
  return(list(Mean = mean,MAE = mae,RMSE = rmse))
}

table_error_3 <- function(benchmark = benchmark_error, model = model_error, ns = dns_results){
  library(xtable)
  maturity <- c(3,6,9,seq(12,360,by=12))
  maturity2 <- c()
  for (i in 1:16){maturity2[i] <- maturity[2*i]}
  table <- 100*cbind(Mean = as.numeric(model$Mean),
                     Mean = as.numeric(benchmark$Mean),Mean = as.numeric(ns$Mean),
                     RMSE = as.numeric(model$RMSE),
                     RMSE = as.numeric(benchmark$RMSE),RMSE = as.numeric(ns$RMSE), MAE = as.numeric(model$MAE), MAE = as.numeric(benchmark$MAE),
                     MAE = as.numeric(ns$MAE),
                     rMAE = as.numeric(model$MAE)/(100*as.numeric(benchmark$MAE)),rMAE2 = as.numeric(model$MAE)/(100*as.numeric(ns$MAE)))
  
  rownames(table) <- maturity
  table2 <- matrix(NA,nr=16,nc=11)
  for (i in 1:16){table2[i,] <-  table[2*i,]}
  rownames(table2) <- maturity2
  #maturity2 <- c(maturity[18:34])
  #table2 <- cbind(table[1:17,],maturity2, table[c(18:33,NA), ])
  m <- c(mean(model$Mean),mean(benchmark$Mean),mean(ns$Mean),mean(model$RMSE),mean(benchmark$RMSE),
         mean(ns$RMSE),mean(model$MAE),mean(benchmark$MAE),mean(ns$MAE),mean(as.numeric(model$MAE)/(100*as.numeric(benchmark$MAE))),
         mean(as.numeric(model$MAE)/(100*as.numeric(ns$MAE))))
  xtable(rbind(table2,100*m))
  
}


# DNS functions stady state-----------------------------------------------------------

self_spec_list <- list()
# We want to specify the H matrix ourselves
self_spec_list$H_spec <- TRUE
# We have got 6 state parameters: 3 factors and 3 fixed means
self_spec_list$state_num <- 6
# In total we need 20 parameters:
# 1 for lambda
# 1 for sigma2 (H)
# 6 for the variance - covariance matrix Sigma_eta (Q)
# 9 for the vector autoregressive coefficient matrix Phi
# 3 for the means mu
self_spec_list$param_num <- 20
# R is a fixed diagonal matrix
self_spec_list$R <- diag(1, 6, 6)
# P_inf is a matrix of zeroes, as all state parameters are stationary
self_spec_list$P_inf <- matrix(0, 6, 6)
# Needed because we want to use collapse = TRUE
# The fixed means only appear in the state equations,
# not in the observation equations. So the 4th, 5th, and 6th state parameters
# are state_only.

self_spec_list$state_only <- 4:6
self_spec_list$sys_mat_fun <- function(param){
  # Maturities of the interest rates
  #maturity <- c(3, 6, 12, 24, 36, 60, 84, 120)
  maturity <- c(3,6,9,seq(12,360,by=12))
  # The constant lambda
  lambda <- exp(2 * param[1])
  # The variance of the observation errors
  sigma2 <- exp(2 * param[2])
  H <- sigma2 * diag(1, NROW(maturity), NROW(maturity))
  # Z matrix corresponding to the factors
  lambda_maturity <- lambda * maturity
  z <- exp(-lambda_maturity)
  Z <- matrix(1, NROW(maturity), 3)
  Z[, 2] <- (1 - z) / lambda_maturity
  #If we want to obtain standard errors of certain (functions of) parameters, we specify transform_fun.
  #Getting proper initial values can be a bit tricky, see vignette("dictionary", "statespacer") for details. Playing
  #around with the parameters corresponding to , , the diagonal elements of , and the diagonal elements
  #of lead to the following initial values:
  #  However, for building this vignette, we make use of the already optimal parameters, to reduce the building
  #time. But you can go ahead and use the initial values as specified above!
  Z[, 3] <- Z[, 2] - z
  # Variance of the state disturbances
  Q <- statespacer::Cholesky(param = param[3:8], decompositions = FALSE, format = matrix(1, 3, 3))
  # Vector autoregressive coefficient matrix, enforcing stationarity
  Tmat <- CoeffARMA(A = array(param[9:17], dim = c(3, 3, 1)),
                    variance = Q,
                    ar = 1, ma = 0)$ar[,,1]
  # Initial uncertainty of the factors
  T_kronecker <- kronecker(Tmat, Tmat)
  Tinv <- solve(diag(1, dim(T_kronecker)[1], dim(T_kronecker)[2]) - T_kronecker)
  vecQ <- matrix(Q)
  vecPstar <- Tinv %*% vecQ
  P_star <- matrix(vecPstar, dim(Tmat)[1], dim(Tmat)[2])
  # Adding parts corresponding to the fixed means to the system matrices
  Z <- cbind(Z, matrix(0, NROW(maturity), 3)) # Not used in the observation equation
  Q <- BlockMatrix(Q, matrix(0, 3, 3)) # Fixed, so no variance in its errors
  a1 <- matrix(param[18:20], 6, 1) # Fixed means go into the initial guess
  Tmat <- cbind(Tmat, diag(1, 3, 3) - Tmat)
  Tmat <- rbind(Tmat, cbind(matrix(0, 3, 3), diag(1, 3, 3)))
  P_star <- BlockMatrix(P_star, matrix(0, 3, 3))
  # Return the system matrices
  return(list(H = H, Z = Z, Tmat = Tmat, Q = Q, a1 = a1, P_star = P_star))
}

self_spec_list$transform_fun <- function(param){
  lambda <- exp(2 * param[1])
  sigma2 <- exp(2 * param[2])
  means <- param[18:20]
  return(c(lambda, sigma2, means))
}
construct_model_parameter_dns <- function(param = dns_param, maturity) 
{
  A      <- param$autoregressiveT   # AR matrix
  MU     <- as.numeric(param$parameter[3:5,2])         # mean vector
  Q      <- param$Variance_of_the_stateQ   # cov in state eq.
  lambda <- as.numeric(param$parameter[1,2])          # lambda
  sigma2 <- as.numeric(param$parameter[2,2])
  H <- sigma2 * diag(1, NROW(maturity), NROW(maturity))
  
  gnk <- 3
  Phi0  <- (diag(gnk)-A)%*%MU 
  delta <- 1/12
  # eval(parameter_load)
  a0 <- param$a0
  P0 <- diag(param$P0)
  npara <- 20
  v.mat <- maturity
  ct <- rep(0,NROW(maturity))
  Zt <- NS.B(lambda,v.mat)
  #h.tr <- sqrt(1/h - 16)
  GGt <- H
  HHt <- Q
  Tt <- A
  dt <- Phi0 
  delta <- 1/12
  
  return(list(ct = ct, Zt = Zt, GGt = GGt, dt = dt, Tt = Tt, 
              HHt = HHt, maturity = maturity, a0 = a0, P0 = P0, delta_Temps = delta))
}

construct_model_parameter_dns_logA0 <- function(param = dns_param, maturity) {
  
  p = 3
  delta <- 1/12
  nb_Maturity <- length(maturity)
  eval(parameter_load)
  #state_mean <- theta_p
  #state_var <- sigma^2 %*% solve(1 - (1 - kappa_p)^2)
  #a0 <- state_mean
  #P0 <- state_var
  a0 <- param$a0
  P0 <- diag(param$P0)
  #lambda <- lambda* delta_Temps # add this line for test
  logAA <- 0 #test DNS
  
  BB_1 <- delta * maturity
  BB_2 <- delta*(1-exp(-lambda*maturity))/lambda    # test  DNS
  BB_3 <- delta*(1-exp(-lambda*maturity))/lambda  - delta*maturity*exp(-lambda*maturity)#test DNS
  aa <- -logAA/(delta*maturity)
  BB <- matrix(c(BB_1/(delta*maturity), BB_2/(delta*maturity), BB_3/(delta*maturity)),nc = 3, byrow = F)
  ct <- aa
  Zt <- BB
  #h.tr <- sqrt(1/h - 16)
  GGt <- diag(rep(h, nb_Maturity))
  HHt <- t(sigma) %*% rho_fin %*% sigma
  Tt <- diag(3) - kappa_p
  #dt <- kappa_p %*% theta_p
  dt <- kappa_q %*% theta_q
  
  
  return(list(ct = ct, Zt = Zt, GGt = GGt, dt = dt, Tt = Tt, 
              HHt = HHt, a0 = a0, P0 = P0, maturity = maturity, delta_Temps = delta))
}

g3d_log_like_dns <- function (par=dns_param, y, maturity, model_type, delta_Temps, constraints, 
                              a0 = NULL, P0 = NULL, output = "loglik", Type = "discrete",par_input="real") 
{
  
  if(par_input=="working"){
    param$autoregressiveT <- matrix(par[1:9],nc=3); param$Variance_of_the_stateQ <- matrix(c(par[10:12],par[11],par[13:14],par[12],par[14:15]),nc=3);
    param$parameter = matrix(c("lambada","sigma2","mu1","mu2","mu3",par[16:20]),nr=5, nc = 2)
    par <- param
  }
  p=3
  n <- ncol(y)
  d <- nrow(y)
  
  
  liste_Param_Espace_Etat <- construct_model_parameter_dns(par, maturity = maturity)
  
  if (is.null(a0)){ 
    a0 <- liste_Param_Espace_Etat$a0}else{a0 <- a0[1:p]}
  
  if(is.null(P0)){ 
    P0 <- liste_Param_Espace_Etat$P0}else{P0 <- as.matrix(as.matrix(P0)[1:p, 1:p])}
  
  dt <- liste_Param_Espace_Etat$dt
  Zt <- liste_Param_Espace_Etat$Zt
  ct <- liste_Param_Espace_Etat$ct
  Tt <- liste_Param_Espace_Etat$Tt
  HHt <- liste_Param_Espace_Etat$HHt
  GGt <- liste_Param_Espace_Etat$GGt
  eigen_Values <- try(eigen(HHt, only.values = TRUE), silent = TRUE)
  #if (is.null(attr(eigen_Values, "class")) == FALSE) {
  #  return(10^10)
  #}else {
  # bool_Eigen_Neg <- any(eigen_Values$values < 0)
  #  if (bool_Eigen_Neg == TRUE) 
  #   return(10^10)
  #}
  
  
  out <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, 
             Zt = Zt, HHt = HHt, GGt = GGt, yt = y)
  if(output == "loglik"){
    logLik <- out$logLik
    if(is.na(logLik) == TRUE){
      return(10^10)
    }else{return(-logLik)} 
  }else if(output == "full"){
    states_filter_mean <- out$att
    states_filter_var <- out$Ptt
    states_pred_mean <- out$at
    states_pred_var <- out$Pt
    states_intercept_a <- out$ct
    states_factor_B <- out$Zt
    error_pred_mean <- out$vt
    states_smooth_mean <- states_filter_mean
    states_smooth_var <- states_filter_var
    for (i in ncol(y):2) {
      J <- states_filter_var[, , i - 1] %*% t(Tt) %*% 
        solve(states_pred_var[, , i])
      states_smooth_mean[, i - 1] <- states_filter_mean[, 
                                                        i - 1] + J %*% (states_smooth_mean[, i] - states_pred_mean[, 
                                                                                                                   i])
      states_smooth_var[, , i - 1] <- states_filter_var[, 
                                                        , i - 1] + J %*% (states_smooth_var[, , i] - 
                                                                            states_pred_var[, , i]) %*% t(J)
    }
    
    states_filter_mean <- data.frame(t(states_filter_mean))
    dimnames(states_filter_mean)[[2]] <- paste0("x", 1:p, 
                                                "_filter")
    
    states_filter_var <- data.frame(t(apply(states_filter_var, 
                                            3, diag)))
    
    
    dimnames(states_filter_var)[[2]] <- paste0("x", 1:p, 
                                               "_filter_var")
    states_pred_mean <- data.frame(t(states_pred_mean))
    dimnames(states_pred_mean)[[2]] <- paste0("x", 1:p, 
                                              "_pred")
    
    states_pred_var <- data.frame(t(apply(states_pred_var, 
                                          3, diag)))
    
    
    dimnames(states_pred_var)[[2]] <- paste0("x", 1:p, "_pred_var")
    states_smooth_mean <- data.frame(t(states_smooth_mean))
    dimnames(states_smooth_mean)[[2]] <- paste0("x", 1:p, 
                                                "_smooth")
    
    states_smooth_var <- data.frame(t(apply(states_smooth_var, 
                                            3, diag)))
    
    
    dimnames(states_smooth_var)[[2]] <- paste0("x", 1:p, 
                                               "_smooth_var")
    ticker <- c(date_Milieu_Mois, NA)
    states_filter_mean <- rbind(states_filter_mean, NA)
    states_filter_var <- rbind(states_filter_var, NA)
    states_smooth_mean <- rbind(states_smooth_mean, NA)
    states_smooth_var <- rbind(states_smooth_var, NA)
    results <- data.frame(ticker, states_filter_mean, states_filter_var, 
                          states_pred_mean, states_pred_var, states_smooth_mean, 
                          states_smooth_var)
    out$summary <- results
    error_pred_mean <- data.frame(ticker[-length(ticker)], 
                                  t(error_pred_mean))
    dimnames(error_pred_mean)[[2]] <- c("date", paste0("maturity_", 
                                                       maturity_all_years))
    out$error_pred <- error_pred_mean
    out$aa <- ct
    out$B <-  states_factor_B
    out$ht <- HHt
    return(out)
  }else{"Wrong Type of Output"}
}

work_param_dns <- function(par){
  param <- list()
  param$autoregressiveT <- matrix(par[1:9],nc=3); param$Variance_of_the_stateQ <- matrix(c(par[10:12],par[11],par[13:14],par[12],par[14:15]),nc=3);
  param$parameter = matrix(c("lambada","sigma2","mu1","mu2","mu3",par[16:20]),nr=5, nc = 2)
  param$a0 = par[21:23]; param$P0 <- par[24:26]
  return(param)
}

NS.B <- function(lambda, tau)
{
  col1 <- rep.int(1,length(tau))
  col2 <- (1-exp(-lambda*tau))/(lambda*tau)
  col3 <- col2-exp(-lambda*tau) 
  return(cbind(col1,col2,col3))
}

creat_output_file_pred <- function(type){
  require(tidyverse)
  date <- c("2017-01-01","2018-01-01","2019-01-01","2020-01-01","2021-01-01","2022-02-10")
  start.date <- as.Date("1986-01-01")
  end.date <- as.Date("2022-01-31")
  model_type <- "discrete"
  constraints <- list(lambda1=0, 
                      lambda2="include", 
                      rho="include",
                      p=3)
  eval(load_data2)
  data.set <- cbind.data.frame(date_Milieu_Mois,t(y))
  load(paste0(getwd(),"/all_three_models_parameters.RData")) # loading parameters
  load(paste0(getwd(),"/dafns_test_dns.RData"))
  load(paste0(getwd(),"/DNS_out_of_sample_parameters.RData"))
  all_models$dns <- dafns_test
  Y <- data.set %>% dplyr::filter(date[1] < data.set[,1])
  predict <- list()
  if(type == "DTAFNS"){
    for(j in 1:(length(date)-1)){
      param <- list()
      param[[j]] <- work_param(unlist(all_models$dafns[[2]][j]))
      #param[[j]] <- work_param(unlist(parameter[[j]]))
      train <- data.set %>% dplyr::filter(data.set[,1] < date[j])
      validation <- data.set %>% dplyr::filter(date[j] < data.set[,1] & data.set[,1] <= date[j+1])
      y1 <- train[,-1]
      y2 <- validation[,-1]
      date_Milieu_Mois <- train[,1]
      results <- g3d_log_like_2(param[[j]], y=t(y1), maturity, model_type, delta_Temps, constraints,
                                par_input="real", output="full")
      a0 <- as.numeric(results$summary[NROW(results$summary)-1,2:4])
      P0 <- diag(as.numeric(results$summary[NROW(results$summary)-1,5:7]))
      date_Milieu_Mois <- train[,1]
      liste_Param_Espace_Etat <- construct_model_parameter(param =  param[[j]], maturity = maturity)
      dt <- liste_Param_Espace_Etat$dt
      Zt <- liste_Param_Espace_Etat$Zt
      ct <- liste_Param_Espace_Etat$ct
      Tt <- liste_Param_Espace_Etat$Tt
      HHt <- liste_Param_Espace_Etat$HHt
      GGt <- liste_Param_Espace_Etat$GGt
      
      date_Milieu_Mois <- validation[,1]
      results <- g3d_log_like_2(param[[j]], t(y2), maturity, model_type, delta_Temps, constraints,
                                a0=a0, P0 = P0,
                                par_input="real", output="full")
      yfitted <- matrix(nrow=nrow(t(y2)), ncol=ncol(t(y2)))
      dffact <- results$summary
      ind <- which(is.na(dffact$ticker)==FALSE) #index of dates that are not NA
      
      DateFactorsG3 <- as.Date(as.character(dffact$ticker[ind]))
      SmoothedFactorsG3 <- cbind(as.numeric(as.character(dffact$x1_smooth[ind])), 
                                 as.numeric(as.character(dffact$x2_smooth[ind])),
                                 as.numeric(as.character(dffact$x3_smooth[ind])))
      rownames(SmoothedFactorsG3) <- as.character(DateFactorsG3)
      
      for(i in 1:ncol(yfitted)) yfitted[,i] <- ct + Zt %*% SmoothedFactorsG3[i,]
      predict$dtafns[[j]] <- yfitted
    }
  }else if(type == "DG3"){
    for(j in 1:(length(date)-1)){
      #a0 <- c(0.04, 0.03, 0.03) 
      #P0 <- diag(c(0.002^2, 0.002^2, 0.002^2))
      param <- list()
      param[[j]] <- unlist(all_models$dg3[[2]][j])
      train <- data.set %>% dplyr::filter(data.set[,1] < date[j])
      validation <- data.set %>% dplyr::filter(date[j] < data.set[,1] & data.set[,1] <= date[j+1])
      y1 <- train[,-1]
      y2 <- validation[,-1]
      date_Milieu_Mois <- train[,1]
      results <- g3d_Log_Lik(param[[j]], t(y1), maturity, model_type, delta_Temps, constraints,
                             par_input="working", output="full")
      a0 <- as.numeric(results$summary[NROW(results$summary)-1,2:4])
      P0 <- diag(as.numeric(results$summary[NROW(results$summary)-1,5:7]))
      param2 <- working_to_real(param[[j]], model_type, delta_Temps, constraints)
      liste_Param_Espace_Etat <- construct_Model_g3d(param =  param2, maturity = maturity)
      dt <- liste_Param_Espace_Etat$dt
      Zt <- liste_Param_Espace_Etat$Zt
      ct <- liste_Param_Espace_Etat$ct
      Tt <- liste_Param_Espace_Etat$Tt
      HHt <- liste_Param_Espace_Etat$HHt
      GGt <- liste_Param_Espace_Etat$GGt
      
      date_Milieu_Mois <- validation[,1]
      results <- g3d_Log_Lik(param[[j]], t(y2), maturity, model_type, delta_Temps, constraints,
                             par_input="working", output="full",
                             a0=a0,P0=P0)
      #a0 <- as.numeric(results$summary[NROW(results$summary)-1,2:4])
      yfitted <- matrix(nrow=nrow(t(y2)), ncol=ncol(t(y2)))
      dffact <- results$summary
      ind <- which(is.na(dffact$ticker)==FALSE) #index of dates that are not NA
      
      DateFactorsG3 <- as.Date(as.character(dffact$ticker[ind]))
      SmoothedFactorsG3 <- cbind(as.numeric(as.character(dffact$x1_smooth[ind])), 
                                 as.numeric(as.character(dffact$x2_smooth[ind])),
                                 as.numeric(as.character(dffact$x3_smooth[ind])))
      rownames(SmoothedFactorsG3) <- as.character(DateFactorsG3)
      
      for(i in 1:ncol(yfitted)) yfitted[,i] <- ct + Zt %*% SmoothedFactorsG3[i,]
      predict$dg3[[j]] <- yfitted
    }
  }else if(type == "DNS"){
    for(j in 1:(length(date)-1)){
      param <- list()
      param[[j]] <- work_param(unlist(parameter[j]))   #New DNS parameters for logAt=0
      train <- data.set %>% dplyr::filter(data.set[,1] < date[j])
      validation <- data.set %>% dplyr::filter(date[j] < data.set[,1] & data.set[,1] <= date[j+1])
      y1 <- train[,-1]
      y2 <- validation[,-1]
      date_Milieu_Mois <- train[,1]
      results <- g3d_log_like_2(param[[j]], y=t(y1), maturity, model_type, delta_Temps, constraints,
                                par_input="real", output="full")
      a0 <- as.numeric(results$summary[NROW(results$summary)-1,2:4])
      P0 <- diag(as.numeric(results$summary[NROW(results$summary)-1,5:7]))
      date_Milieu_Mois <- train[,1]
      liste_Param_Espace_Etat <- construct_model_parameter_dns_logA0(param =  param[[j]], maturity = maturity)
      dt <- liste_Param_Espace_Etat$dt
      Zt <- liste_Param_Espace_Etat$Zt
      ct <- liste_Param_Espace_Etat$ct
      Tt <- liste_Param_Espace_Etat$Tt
      HHt <- liste_Param_Espace_Etat$HHt
      GGt <- liste_Param_Espace_Etat$GGt
      
      date_Milieu_Mois <- validation[,1]
      # for results construct_model_parameter should check to be sure about logA = o and BB
      results <- g3d_log_like_2(param[[j]], t(y2), maturity, model_type, delta_Temps, constraints,
                                a0=a0, P0 = P0,
                                par_input="real", output="full")
      yfitted <- matrix(nrow=nrow(t(y2)), ncol=ncol(t(y2)))
      dffact <- results$summary
      ind <- which(is.na(dffact$ticker)==FALSE) #index of dates that are not NA
      
      DateFactorsG3 <- as.Date(as.character(dffact$ticker[ind]))
      SmoothedFactorsG3 <- cbind(as.numeric(as.character(dffact$x1_smooth[ind])), 
                                 as.numeric(as.character(dffact$x2_smooth[ind])),
                                 as.numeric(as.character(dffact$x3_smooth[ind])))
      rownames(SmoothedFactorsG3) <- as.character(DateFactorsG3)
      
      for(i in 1:ncol(yfitted)) yfitted[,i] <- ct + Zt %*% SmoothedFactorsG3[i,]
      predict$dns[[j]] <- yfitted
    }
  }else{print("wrong input")}
  
  return(list(predict, Y = Y))
}
