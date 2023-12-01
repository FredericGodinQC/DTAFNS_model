



# Codes for the second paper --------------------------------------------------

source(paste0(getwd(),"/function_DTAFNS.R"), echo=TRUE)

#Step 3: Set path
path <- 'D:/Concordia Courses/Research paper/Research/Second/Codes_myself'
setwd(path)

############################################
############################################
#### Model, data and optimization input ####
############################################
############################################

#############
#Model input#
#############

delta_Temps <- 1/12


#############
#Data input #
#############
start.date <- as.Date("1986-01-01")  # "1986-01-01"
end.date <- as.Date("2022-01-31")  # "2018-10-31"  "2022-02-09"

maturity_type <- "all_years" #"all" or "all_years" or "BoC"
#"all": all ZCB maturities in the data set are used -> gives NA (seems like a bug in fkf)
#"all_years": all integer maturities are included + short-end curve (0.25, 0.50, 0.75)
#"BoC": ZCB maturities used in the Bank of Canada working paper (Working Paper 2001-15)
#Bank of Canada working paper uses the following maturities:
#maturity_BoC <- c(0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30)

#############################
#############################
#### Load data           ####
#############################
#############################


eval(load_data2)



######################################
######################################
#### Parameter estimation manual  ####
######################################
######################################

#initial values

param <- list(theta_q = c(0.2, -0.1, 0.22), sigma=c(0.00524, 0.00482, 0.00788),
              lambda = 0.04,
              rho = c(0.140, -0.786, -0.542),
              eta = c(0, 0.001, 10), h = 0.003, 
              a0 = c(0.05, 0.04, 0.03),
              P0 = c(0.002^2, 0.002^2, 0.002^2),
              delta_Temps=1/12)

a0 <- c(0.05, 0.04, 0.03) 
P0 <- diag(c(0.002^2, 0.002^2, 0.002^2))


par <- param_work(param)
#param$eta[1] <- 0
#par <- Optim$pars

#log-likelihood at initial values
loglik0 <- -g3d_log_like_2(par, y, maturity,
                           delta_Temps, constraints, #a0=a0, P0=P0,
                           par_input="working",output="loglik", Type = "discrete")


#optimization

ineqlb1 <- c(rep(-1,length(maturity)),rep(-1,length(maturity)),
             rep(-1,length(maturity)))
inequb1 <-  c(rep(1,length(maturity)),rep(1,length(maturity)),
              rep(1,length(maturity)))
start.time <- proc.time()["elapsed"]


Optim <- solnp(par, fun=g3d_log_like_2,y = y,
               maturity=maturity,
               delta_Temps=delta_Temps, constraints=constraints,
               a0=a0, P0=P0,
               output="loglik",par_input="working",
               control=list(trace=0, tol=1e-8),
               LB = c(rep(-1000,3),rep(0,4),rep(-1,3),rep(-1000,3),0,rep(-1000,3),rep(0,4)),
               UB = c(rep(1000,3),rep(1000,3),100,rep(1,3),rep(1000,3),1000,rep(1000,3),rep(1000,3),1),
               ineqfun = ineqn1, ineqLB = ineqlb1,
               ineqUB = inequb1
)



# Test for DNS
save(Optim, file = paste0(getwd(),"/DNS_test_Optim_solnp.RData"))
load(paste0(getwd(),"/Optim_solnp_new.RData"))
#load(paste0(getwd(),"/Optim_solnp_2018.RData"))
#param <- work_param(Optim_go$pars)

end.time <- proc.time()["elapsed"]
run.time <- as.vector(end.time - start.time)

#results
param <- work_param(Optim$pars)
param$theta_q[1] <- 0
loglik <- -Optim$values[length(Optim$values)]
loglik; param
results <- g3d_log_like_2(param, y, maturity, delta_Temps, constraints,
                          #a0=a0, P0=P0,
                          par_input="real", output="full")


filename <- paste(start.date, "_",
                  "G", constraints$p, "_", model_type, 
                  "_lambda1_", constraints$lambda1,
                  "_lambda2_", constraints$lambda2,
                  "_rho_", constraints$rho,
                  sep="")

write.csv(results$summary, paste0(filename, "_summary", ".csv"))
write.csv(results$error_pred, paste0(filename, "_error_pred", ".csv"))

#####################################################
#Figures to evaluate model fit                      #
#####################################################

#####Bank of Canada yields#####
# date_Milieu_Mois
# maturity
# y

#convert dates into a numeric vector (useful for figures)
dates <- date_Milieu_Mois
t_ <- as.POSIXlt(dates, format="%Y-%m-%d")
n.days <- rep(365, length(t_))
n.days[((1900 + t_$year) %% 4)==0] <- 366 #leap years
t_  <- 1900 + t_$year + (t_$yday+1)/n.days
#when yday=0, we are jan. 1st
#as.POSIXlt("2012-01-01",format="%Y-%m-%d")$yday
#####G3 model#########
dffact <- results$summary
ind <- which(is.na(dffact$ticker)==FALSE) #index of dates that are not NA

DateFactorsG3 <- as.Date(as.character(dffact$ticker[ind]))
SmoothedFactorsG3 <- cbind(as.numeric(as.character(dffact$x1_smooth[ind])), 
                           as.numeric(as.character(dffact$x2_smooth[ind])),
                           as.numeric(as.character(dffact$x3_smooth[ind])))
rownames(SmoothedFactorsG3) <- as.character(DateFactorsG3)

#short rate
ShortRateG3 <- rowSums(SmoothedFactorsG3) - SmoothedFactorsG3[, 3] 
names(ShortRateG3) <- DateFactorsG3

#fitted yields
#########################################################################
#### Doc 1                                                           ####
#### state space system in package fkf                               ####
#### alpha_t+1 = d_t + T_t * alpha_t + eta_t, eta ~ N(0, HHt)        ####
#### y_t = c_t + Z_t * alpha_t + e_t, e_t ~ N(0, GGt)                ####
#### state space system in paper notation                            ####
#### x_t+1 = b + Dx_t + z_t+1, z_t+1 ~ N(0,Q)                        ####
#### y_t = a + B * x_t + eta_t, eta_t ~ N(0, H)                      ####
#########################################################################
liste_Param_Espace_Etat <- construct_model_parameter(param = param, maturity = maturity)
dt <- liste_Param_Espace_Etat$dt
Zt <- liste_Param_Espace_Etat$Zt
ct <- liste_Param_Espace_Etat$ct
Tt <- liste_Param_Espace_Etat$Tt
HHt <- liste_Param_Espace_Etat$HHt
GGt <- liste_Param_Espace_Etat$GGt

yfitted <- matrix(nrow=nrow(y), ncol=ncol(y))


for(i in 1:ncol(yfitted)) yfitted[,i] <- ct + Zt %*% SmoothedFactorsG3[i,]
mymodel_last_year <- yfitted[,ncol(yfitted)]

#--------------------
# figure parameters
#--------------------

xlim <- range(t_)
t.at <- seq(from=1990,to=2021,by=5)

cex.mtext  <- 0.85
cex.axis   <- 0.85
cex.legend <- 0.70

ylab.line <- 0.5
height <- 2.5

###########################################################################
#FIGURE 1: Model-implied factors and short rate                           #
###########################################################################

filename <- "Figure1"
pdf(file=paste(filename,"pdf",sep="."), width=6.5, height=height)
par(mfrow=c(1,1), mar=c(2.5, 2.75, 0.5, 0.5))

ylim.min <- min(SmoothedFactorsG3, ShortRateG3)
ylim.max <- max(SmoothedFactorsG3, ShortRateG3)
ylim <- c(ylim.min, ylim.max)

legend.name <- c("Factor 1", "Factor 2", "Factor 3", "Short rate")
col.all <- c("blue", "green4", "red", "black")

lty.all <- c(5,4,3,1)
lwd.all <- c(2,2,3,2)
pch.all <- rep(NA,length(legend.name))
pt.cex.all <- rep(NA,length(legend.name))
legend.ncol <- length(legend.name)

plot(t_, ShortRateG3, type="n", main="", las=1,
     xaxs="i", xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=xlim, ylim=ylim)#, log="y")

for(i in 1:3) {
  lines(t_, SmoothedFactorsG3[,i], col=col.all[i], lty=lty.all[i], lwd=lwd.all[i])
}
lines(t_, ShortRateG3, col=col.all[4], lty=lty.all[4], lwd=lwd.all[4])
axis(side=1, at=t.at, cex.axis=cex.axis, padj=-0.6)
axis(side=2, cex.axis=cex.axis, las=1)
legend('topright', inset=c(0.02,0), legend=legend.name, lty=lty.all, lwd=lwd.all,
       col=col.all, pch=pch.all, pt.cex=pt.cex.all, bty="n", cex=cex.legend, #horiz=TRUE,
       ncol=legend.ncol, x.intersp=.5, y.intersp=1.5, seg.len=2.8)
dev.off()

###########################################################################
#FIGURE B.1: Model-implied versus observed rates                          #
###########################################################################

filename <- "FigureB1"
pdf(file=paste(filename,"pdf",sep="."), width=6.5, height=height)
par(mfrow=c(1,1), mar=c(2.5, 2.75, 0.5, 0.5))

rates <- c(3, 10*12)
yrates <- t(y[which(maturity %in% rates),])
yfittedrates <- t(yfitted[which(maturity %in% rates),])
yrates[is.na(yrates)==T] <- 0   # man ezafe kardam
ylim.min <- min(yrates, yfittedrates)
ylim.max <- max(yrates, yfittedrates)
ylim <- c(ylim.min, ylim.max)

legend.name <- c("3-month model", "3-month obs", "10-year model", "10-year obs")
col.all <- c("black", "black", "blue", "blue")

lty.all <- c(1,NA,5,NA)
lwd.all <- c(1,NA,1,NA)
pch.all <- c(NA,1,NA,2)
pt.cex.all <- c(NA,0.4,NA,0.3)
legend.ncol <- length(legend.name)

plot(t_, yrates[,1], type="n", main="", las=1,
     xaxs="i", xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=xlim, ylim=ylim)#, log="y")

abline(h=0, lty=2, col="black", lwd=0.5)
for(i in 1:ncol(yrates)) {
  lines(t_, yfittedrates[,i], col=col.all[2*i-1], lty=lty.all[2*i-1], lwd=lwd.all[2*i-1])
  points(t_, yrates[,i], col=col.all[2*i], pch=pch.all[2*i], cex=pt.cex.all[2*i])
}
axis(side=1, at=t.at, cex.axis=cex.axis, padj=-0.6)
axis(side=2, at=c(0,0.05,0.1), cex.axis=cex.axis, las=1)
legend('topright', inset=c(0.02,0), legend=legend.name, lty=lty.all, lwd=lwd.all,
       col=col.all, pch=pch.all, pt.cex=pt.cex.all, bty="n", cex=cex.legend, #horiz=TRUE,
       ncol=legend.ncol, x.intersp=.5, y.intersp=1.5, seg.len=2)

dev.off()

###########################################################################
#FIGURE B.2: Model-implied and observed yield curves                      #
###########################################################################

dates_select <- c("2006-12-29", "2008-12-31", "2016-06-30", "2018-10-31")
dates_i <- which(as.character(dates) %in% dates_select)
graph_per_row <- 2 #number of graphs per row

#--------------------
# figure parameters
#--------------------
cex.mtext  <- 0.75
cex.axis   <- 0.85
cex.legend <- 0.70
ylab.line <- 0.5

height <- 2.5*ceiling(length(dates_i)/graph_per_row)

filename <- "FigureB3"
pdf(file=paste(filename,"pdf",sep="."), width=6.5, height=height)
par(mfrow=c(ceiling(length(dates_i)/graph_per_row), graph_per_row), 
    mar=c(2.25, 2.75, 1.75, 0.5))

for(i in dates_i){
  
  x <- maturity/12
  xlim <- c(0,30)
  
  lty.all <- c(1,NA)
  lwd.all <- c(1,NA)
  pch.all <- c(NA,1)
  pt.cex.all <- c(NA,0.4)
  
  legend.name <- c("Model-implied curve", "Observed curve")
  legend.ncol <- length(legend.name)
  
  ylim.min <- 0     #min(y[,i], yfitted[,i], na.rm=TRUE)
  ylim.max <- 0.06
  #ylim.max <- 0.042 #0.135 #max(y[,i], yfitted[,i], na.rm=TRUE)
  ylim <- c(ylim.min, ylim.max)
  
  plot(x, yfitted[,i], type="n", main="", las=1,
       xaxs="i", xlab="", ylab="", xaxt="n", yaxt="n",
       xlim=xlim, ylim=ylim)
  lines(x, yfitted[,i], col="green", lty=lty.all[1], lwd=lwd.all[1])
  points(x, y[,i], col="black", pch=pch.all[2], cex=pt.cex.all[2],lwd=3)
  #lines(x, yfitted2[,i], col="blue", lty=2, lwd=2)
  #lines(x, yfitted_dns[,i], col="red", lty=4, lwd=3)
  axis(side=1, at=seq(from=0,to=30,by=5), cex.axis=cex.axis, padj=-0.6)
  axis(side=2, at=c(0,0.02,0.04,0.06), cex.axis=cex.axis, las=1) # 0.05
  mtext(dates[i], side=3, line=ylab.line, las=0, cex=cex.mtext)
  legend("topleft",legend = c("Observed yield curve","DTAFNS","DG3", "DNS"),
         bg = rgb(0, 0, 0, alpha = 0.02),
         cex = 0.6 ,col=c("black","green","blue","red"),
         lty = c(3,1, 2, 4))
  
}

dev.off()
library(xtable)
param_p <- convert_q_to_p(param)

param_total <- rbind(param_p$rho,param_p$kappa_p,param_p$sigma,param_p$eta)
xtable(param_total, digits = 5)



