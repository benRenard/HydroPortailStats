library(HydroPortailStats)

# Select case
folder='HBay_Config'
case='Nyons'

# Import config files and run R-HBay
config=Import_HBayConfig(file.path(folder,case))

H3=Hydro3_HBay(y=config$y,dist=config$dist,prior=config$prior,
               SystErrorIndex=config$SystErrorIndex,
               SystErrorPrior=config$SystErrorPrior,
               options=config$options,
               mcmcoptions=config$mcmcoptions)
Hydro3_Plot(H3)

# Compare with Exe-HBay
D=read.table(file.path(folder,case,'quantiles.bay_qtl'),header = TRUE)

plot(H3$quantile$T,H3$quantile$IC.high,type='l',log='x')
lines(H3$quantile$T,H3$quantile$IC.low)
lines(H3$quantile$T,H3$quantile$q,lwd=2)

lines(D$T,D$IC_lower,lty=2,col='red')
lines(D$T,D$IC_upper,lty=2,col='red')
lines(D$T,D$modal_Q.p.,lty=2,col='red',lwd=2)
