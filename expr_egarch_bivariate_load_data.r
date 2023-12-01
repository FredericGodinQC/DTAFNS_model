#-------------------------------------------------------------------------------
#LOAD DATA FOR EGARCH BIVARIATE
#-------------------------------------------------------------------------------

expr_load_data_egarch <- expression({
                            
############EQUITY INDEX RETURNS#############

DATA_returns <- read.csv("DATA_returns_all.csv")

#TSX return
logreturnTSX <- DATA_returns[["TSX"]]
names(logreturnTSX) <- DATA_returns[["date"]]

#SP500 return
logreturnSP500 <- DATA_returns[["SP500"]]
names(logreturnSP500) <- DATA_returns[["date"]]

logR <- cbind(TSX=logreturnTSX, SP500=logreturnSP500)

############SHORT RATE#############

#G3 model
dffact <- read.csv2("1986-01-01_G3_discrete_lambda1_0_lambda2_include_rho_include_summary.csv", header = TRUE, sep = ",")
dffact <- subset(dffact, select=-c(X))

ind <- which(is.na(dffact$ticker)==FALSE) #index of dates that are not NA

DateFactorsG3 <- as.Date(as.character(dffact$ticker[ind]))
SmoothedFactorsG3 <- cbind(as.numeric(as.character(dffact$x1_smooth[ind])), 
                          as.numeric(as.character(dffact$x2_smooth[ind])),
                          as.numeric(as.character(dffact$x3_smooth[ind])))
rownames(SmoothedFactorsG3) <- as.character(DateFactorsG3)

ShortRateG3 <- rowSums(SmoothedFactorsG3) - SmoothedFactorsG3[,3]
names(ShortRateG3) <- DateFactorsG3

})