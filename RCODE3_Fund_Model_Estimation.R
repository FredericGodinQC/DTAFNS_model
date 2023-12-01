###############################################################################################
#### ESTIMATION OF MIXED BOND AND EQUITY FUND MODEL                                       #####
###############################################################################################

#Install the following packages
library(rugarch)

#Set path 
#path <- 'C:/Dropbox/MaciejEmmanuel/Rcode/WEBSITE/3_Fund_model'
#path <- 'D:/Concordia Courses/Research paper/Research/Second/R_codes_mixed_bond/3_Fund_model'
path <- 'D:/Concordia Courses/Research paper/Research/Second/Codes_myself/3_Fund_model_mymodel'
setwd(path)

###############################################################################################
#####################################  LOAD DATA  #############################################
###############################################################################################

################################
#Load interest rate factor data#
################################

#kQ (under Q) parameters estimated for G3 model
kQ_G3 <- c(0, 0.02332395, 0.02332395) #K_ii

############TERM STRUCTURE SMOOTHED FACTORS#############
dffact <- read.csv2("1986-01-01_G3_discrete_lambda1_0_lambda2_include_rho_include_summary.csv", header = TRUE, sep = ",")
dffact <- subset(dffact, select=-c(X))

ind <- which(is.na(dffact$ticker)==FALSE) #index of dates that are not NA

DateFactorsG3 <- as.Date(as.character(dffact$ticker[ind]))
SmoothedFactorsG3 <- cbind(as.numeric(as.character(dffact$x1_smooth[ind])), 
                           as.numeric(as.character(dffact$x2_smooth[ind])),
                           as.numeric(as.character(dffact$x3_smooth[ind])))
rownames(SmoothedFactorsG3) <- as.character(DateFactorsG3)

#short rate
ShortRateG3 <- rowSums(SmoothedFactorsG3) - SmoothedFactorsG3[, 3]
names(ShortRateG3) <- DateFactorsG3

#factor variations
diffSmoothedFactorsG3 <- diff(SmoothedFactorsG3)                         
#SmoothedFactorsG3[-1,] - SmoothedFactorsG3[-nrow(SmoothedFactorsG3),] - diffSmoothedFactorsG3

#factor variations that incorporate the (1-kQ) term     x_{t+1} - (1-kQ_G3)x_t +((1-kQ_G3[3])^(maturity-1)-1)x^(3)_t
diffSmoothedFactorsG3mod <- SmoothedFactorsG3[-1,] - 
  t((1-kQ_G3)*t(SmoothedFactorsG3[-nrow(SmoothedFactorsG3),])) #+ ((1-kQ_G3[3])^(maturity-1)-1)*SmoothedFactorsG3[-nrow(SmoothedFactorsG3),3]

###################
#load returns data#
###################
DATA_returns <- read.csv("DATA_returns_all.csv")

###########################################################################
#############################   ESTIMATION   ##############################
###########################################################################

######## SELECT FUND #######
#1 = RBC Bond GIF Series 1
#2 = Assumption/Louisbourg Balanced Fund A
#3 = iShares Core Canadian Universe Bond Index ETF (XBB)

fundno = 2

#-----------
#model entry
#-----------
nb_factors <- 3
include_TSX <- TRUE
include_SP500 <- TRUE
volatility <- 'eGARCH'
arOrder <- 0
distribution.model <- "norm"
#valid choices are ?norm? for the normal distibution,
#?snorm? for the skew-normal distribution,
#?std? for the student-t, ?sstd? for the skew-student,
#?ged? for the generalized error distribution, ?sged? for the skew-generalized error distribution,
#?nig? for the normal inverse gaussian distribution,
#?ghyp? for the Generalized Hyperbolic, and
#?jsu? for Johnson's SU distribution.
include.mean <- TRUE
include.1_kQ <- TRUE

#-----------
#load vars
#-----------
armaOrder <- c(arOrder,0)

#load term structure variables
if(nb_factors==3){
  model <- "G3"
  rf_return <- ShortRateG3/12   
  if(include.1_kQ==TRUE){
    diffSmoothedFactors <- cbind(diffSmoothedFactorsG3mod, SmoothedFactorsG3[-nrow(SmoothedFactorsG3),3])   # change diffSmoothedFactorsG3mod
  } else if(include.1_kQ==FALSE){
    diffSmoothedFactors <- diffSmoothedFactorsG3
  }
}

#fund return net of risk-free rate
logreturnfund <- as.numeric(DATA_returns[[3+fundno-1]])    # -1 add
RegYrf <- logreturnfund - rf_return[-length(rf_return)]
names(RegYrf) <- rownames(diffSmoothedFactors)
#remove NA values
ind <- which(is.na(RegYrf)==FALSE)
RegYrf <- RegYrf[ind]
#adjust dimension of term structure variables
diffSmoothedFactors <- diffSmoothedFactors[ind,]

#TSX return
logreturnTSX <- DATA_returns[["TSX"]] - ShortRateG3[-1]/12  
names(logreturnTSX) <- DATA_returns[["date"]]
logreturnTSX <- logreturnTSX[ind]

#SP500 return
logreturnSP500 <- DATA_returns[["SP500"]] - ShortRateG3[-1]/12
names(logreturnSP500) <- DATA_returns[["date"]]
logreturnSP500 <- logreturnSP500[ind]

nb_factors_regression <- nb_factors + 1 # added + 1
  
external.regressors <- NULL
if(nb_factors_regression>0){
  external.regressors <- diffSmoothedFactors[,1:nb_factors_regression]
  if(nb_factors_regression==1) external.regressors <- as.matrix(external.regressors, ncol=1)
} else if(nb_factors_regression==0){
  external.regressors <- NULL   
}
  
if(include_TSX==TRUE){
  external.regressors <- cbind(external.regressors, as.vector(logreturnTSX))
}

if(include_SP500==TRUE){
  external.regressors <- cbind(external.regressors, as.vector(logreturnSP500))
}

if(volatility=="constant") {
  #ARFIMA
  mean.model <- list(armaOrder = armaOrder, 
                     include.mean = include.mean,
                     external.regressors = external.regressors,
                     arfima = FALSE
                     )  
  variance.model <- NULL

  arfima_spec <- arfimaspec( mean.model = mean.model,
                             distribution.model = distribution.model )  

  start.time <- proc.time()["elapsed"]
  bondfundmodel <- try(arfimafit(spec=arfima_spec, data=RegYrf), silent=TRUE)
  #names(bondfundmodel@fit)
  end.time <- proc.time()["elapsed"]
  run.time <- as.vector(end.time - start.time)      
} else {
  #rugarch
  mean.model <- list(armaOrder = armaOrder, 
                     include.mean = include.mean,
                     external.regressors = external.regressors,
                     archm = FALSE,  #Whether to include ARCH volatility in the mean regression
                     archpow = 1,    #Indicates whether to use st.deviation (1) or variance (2) in the ARCH in mean regression
                     arfima = FALSE, #Whether to fractional differencing in the ARMA regression
                     archex = FALSE  #Whether to multiply the last 'archex' external regressors by the conditional standard deviation
                    )  
  variance.model <- list(model = volatility, garchOrder = c(1,1))  
  #variance.model <- list(model = "sGARCH", garchOrder = c(1, 1),
  #        submodel = NULL, external.regressors = NULL, variance.targeting = FALSE)
  #variance.model <- list(model = "eGARCH", garchOrder = c(1, 1),
  #                        external.regressors = as.matrix(logreturnequityindex, ncol=1))
  #variance.model <- list(model = "csGARCH", garchOrder = c(1, 1))
  #variance.model <- list(model = "fGARCH", garchOrder = c(1, 1), submodel="GJRGARCH")
  #variance.model <- list(model = "fGARCH", garchOrder = c(1, 1), submodel="TGARCH")

  rgarch_spec <- ugarchspec( variance.model = variance.model,
                             mean.model = mean.model,
                             distribution.model = distribution.model )
  #fixed.pars=list("alpha1"=0,"beta1"=0)
  #setfixed(rgarch_spec) <- list("alpha1"=0,"beta1"=0)  

  start.time <- proc.time()["elapsed"]
  bondfundmodel <- try(ugarchfit(spec=rgarch_spec, data=RegYrf
                                 #fit.control=list(stationarity=0),
                                 #solver="nloptr",
                                 #solver.control = list(trace = 1)
                                 ), silent=TRUE)
  #names(bondfundmodel@fit)
  end.time <- proc.time()["elapsed"]
  run.time <- as.vector(end.time - start.time)  
}

#results
params <- coef(bondfundmodel)
if("mu" %in% names(params)) theta0 <- as.vector(params["mu"]) else theta0 <- NA
if("ar1" %in% names(params)) ar1 <- as.vector(params["ar1"]) else ar1 <- NA
if("omega" %in% names(params))  omega <- as.vector(params["omega"])  else omega <- NA
if("alpha1" %in% names(params)) a1    <- as.vector(params["alpha1"]) else a1 <- NA
if("beta1" %in% names(params))  b1    <- as.vector(params["beta1"])  else b1 <- NA
if("eta11"  %in% names(params)) g1    <- as.vector(params["eta11"])  else g1 <- NA
if("gamma1" %in% names(params)) g1    <- as.vector(params["gamma1"]) #EGARCH only
if("skew"   %in% names(params)) skew  <- as.vector(params["skew"])   else skew  <- NA
if("shape"  %in% names(params)) shape <- as.vector(params["shape"])  else shape  <- NA
#if("eta21"  %in% names(params)) eta <- as.vector(params["eta21"])  else eta <- NA
#if("lambda" %in% names(params)) pow <- as.vector(params["lambda"]) else pow <- NA

theta1 <- theta2 <- theta3 <- theta_prime <- NA
if(nb_factors_regression==1){
  theta1 <- as.vector(params["mxreg1"])
} else if(nb_factors_regression==2){
  theta1 <- as.vector(params["mxreg1"])
  theta2 <- as.vector(params["mxreg2"])
} else if(nb_factors_regression==3){
  theta1 <- as.vector(params["mxreg1"])
  theta2 <- as.vector(params["mxreg2"])
  theta3 <- as.vector(params["mxreg3"])
} else if(nb_factors_regression==4){           # I add this one
  theta1 <- as.vector(params["mxreg1"])
  theta2 <- as.vector(params["mxreg2"])
  theta3 <- as.vector(params["mxreg3"])
  theta_prime <- as.vector(params["mxreg4"])
}

mxreg <- paste0("mxreg", sum(!is.na(c(theta1,theta2,theta3,theta_prime)))+1)  # I add theta_prime

thetaTSX <- NA
if(include_TSX==TRUE) {
  thetaTSX <- as.vector(params[mxreg])
  mxreg <- paste0("mxreg", sum(!is.na(c(theta1,theta2,theta3,theta_prime)))+2)
}

thetaSP500 <- NA
if(include_SP500==TRUE) thetaSP500 <- as.vector(params[mxreg])

if(volatility=="constant") {
  #parameter in ARFIMA
  sigma <- as.vector(params["sigma"])  
} else {
  sigma <- sqrt(uncvariance(bondfundmodel)) #unconditional stdev of innovation
  #for EGARCH it does: sqrt( exp( as.vector(params["omega"]) / (1-pers) ) )
}


if(length(sigma)==0) sigma <- NA
#sigma <- round(sigma,6)

pers <- bondfundmodel@fit$persistence #persistence in GARCH
if(is.null(pers)) pers <- NA 
pers <- round(pers, 4)
h0_f <- sigma^2 # I added this part for volatility
params <- c(sigma=sigma, theta0=theta0, theta1=theta1, theta2=theta2, theta3=theta3, theta_prime = theta_prime,
            thetaTSX=thetaTSX, thetaSP500=thetaSP500, ar1=ar1,
            skew=skew, shape=shape,  
            omega=omega, a1=a1, b1=b1, g1=g1, pers=pers, h0_f = h0_f)
params_fund <- params
save(params_fund, file = paste0(getwd(),"/fund_param.RData"))
#save(date_Milieu_Mois, file = paste0(getwd(),"/date_Milieu_Mois.RData"))
loglik <- likelihood(bondfundmodel) #bondfundmodel@fit$LLH
if(is.null(loglik)) loglik <- NA
#residuals(bondfundmodel)
#plot(bondfundmodel,which="all")
params_n <- length(coef(bondfundmodel))
if(params_n==0) params_n <- NA
res <- c(fundno=fundno,
         #index=index,
         model=model,  
         factors=nb_factors_regression,
         '1_kQ'=include.1_kQ,
         TSX=include_TSX,
         SP500=include_SP500,
         intercept=include.mean,
         ARorder=arOrder,
         sigmat=volatility,
         dist=distribution.model,
         par=params_n,
         loglik=loglik,
         AIC=2*params_n - 2*loglik,
         BIC=log(length(RegYrf))*params_n - 2*loglik,
         params,
         time=round(run.time,1),
         conv=bondfundmodel@fit$convergence,
         n=length(RegYrf))
as.data.frame(res)

filename <- paste0("___fundno", "_", fundno, "_", "estimation", "_", "single")
write.csv(res, paste(filename,"csv",sep="."))

ticker <- paste0("___fundno", "_", fundno, "_", "diagnostics", "_", "single")
sink(paste0(ticker,".txt"))
show(bondfundmodel)
sink()

