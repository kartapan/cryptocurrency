rm(list=ls()) 
library(dplyr)
library(MSGARCH)
library(rugarch)
library(copula)
library(GAS)
library(RobGARCHBoot)
library(esback)
library(forecast)

setwd("C:/Users/Tapan Kar/Dropbox/Crypto new/Crypto_QF")

data <- read.csv("ethereum data.csv")
close= data$Close
open= data$Open
high= data$High
low= data$Low
l= length(close)

pk_proxy= {(log(high)-log(low))^2}/(4*log(2))
pk_proxy= pk_proxy[-1]

log_ret= c(rep(0,l))

for(i in 2:l){
  log_ret[i]= log(close[i]/close[(i-1)])*100
}

# deleting first observation
log_ret= log_ret[-1]

log_range= c(rep(0,l))

for(i in 1:l){
  log_range[i]= log(high[i])-log(low[(i)])
}

log_range= log_range[-1]

range= sqrt(log_range)



# GARCH Model with sstd
m= 1343
pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

# Garch model forecast

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  fit_mean_garch[q]= as.numeric(pre_mean_garch$pred)
  
  spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sstd")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH)
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_garch_fast1[q]= as.numeric(fore$sigmaFor)
  
  
}

etherum_garch_mean= fit_mean_garch
etherum_garch= pre_garch_fast1

write.csv(etherum_garch_mean,"etherum_garch_mean.csv")
etherum_garch_mean= read.csv("etherum_garch_mean.csv")
etherum_garch_mean= etherum_garch_mean[,-1]
#write.csv(etherum_garch,"etherum_garch.csv")
etherum_garch_vol= read.csv("etherum_garch.csv")
etherum_garch_vol= etherum_garch_vol[,-1]

# NAGARCH Model with sstd

pre_nagarch_binance= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  
  spec.gjrGARCH = ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(1,1),submodel= "NAGARCH"), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sstd")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH,solver = "hybrid")
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_nagarch_binance[q]= as.numeric(fore$sigmaFor)
  
}

bid_nagarch_binance_eth= pre_nagarch_binance
write.csv(bid_nagarch_binance_eth,"bid_nagarch_binance_eth.csv")
nagarch_binance_eth= read.csv("bid_nagarch_binance_eth.csv")
nagarch_binance_eth= nagarch_binance_eth[,-1]


#GJR-GARCH Model with sstd

pre_gjrgarch_fast1=c(rep(0,m))
for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  res_gjrgarch= fit$residuals
  spec.garch = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sstd")
  e_gjrgarch <- ugarchfit(res_gjrgarch, spec=spec.garch)
  e_forecast= ugarchforecast(e_gjrgarch,n.ahead=1,data=res_gjrgarch)
  fore= e_forecast@forecast
  pre_gjrgarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

eth_gjrgarch= pre_gjrgarch_fast1
write.csv(eth_gjrgarch,"eth_gjrgarch.csv")
eth_gjrgarch= read.csv("eth_gjrgarch.csv")
eth_gjrgarch= eth_gjrgarch[,-1]

#  Rolling window forecast for MSGARCH Model
msgarch_var_01_binance= c(rep(0,m))
msgarch_es_01_binance= c(rep(0,m))
msgarch_var_05_binance= c(rep(0,m))
msgarch_es_05_binance= c(rep(0,m))
msgarch_var_025_binance= c(rep(0,m))
msgarch_es_025_binance= c(rep(0,m))
eth_msgarch_vol_binance= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  ms2.garch.n= CreateSpec(variance.spec = list(model = "sGARCH"),
                          distribution.spec = list(distribution = "sstd"),
                          switch.spec = list(K = 2))
  fit.ml <- FitMCMC(spec = ms2.garch.n, data = res_garch)
  ms_garch= predict(fit.ml,nahead = 1)
  eth_msgarch_vol_binance[q]= ms_garch$vol
  
  e_forecast_01= Risk(fit.ml, alpha = 0.01, nahead = 1)
  e_forecast_05= Risk(fit.ml, alpha = 0.05, nahead = 1)
  e_forecast_025= Risk(fit.ml, alpha = 0.025, nahead = 1)
  
  msgarch_var_01_binance[q]= e_forecast_01$VaR
  msgarch_es_01_binance[q]= e_forecast_01$ES
  msgarch_var_05_binance[q]= e_forecast_05$VaR
  msgarch_es_05_binance[q]= e_forecast_05$ES
  
  msgarch_var_025_binance[q]= e_forecast_025$VaR
  msgarch_es_025_binance[q]= e_forecast_025$ES
  
}

eth_msgarch_var_01_binance=  msgarch_var_01_binance
eth_msgarch_es_01_binance=  msgarch_es_01_binance
eth_msgarch_var_05_binance=  msgarch_var_05_binance
eth_msgarch_es_05_binance=  msgarch_es_05_binance
eth_msgarch_vol_binance= eth_msgarch_vol_binance

eth_msgarch_var_025_binance=  msgarch_var_025_binance
eth_msgarch_es_025_binance=  msgarch_es_025_binance

write.csv(eth_msgarch_vol_binance,"eth_msgarch_vol_binance.csv")

write.csv(eth_msgarch_var_01_binance,"eth_msgarch_var_01_binance.csv")
write.csv(eth_msgarch_es_01_binance,"eth_msgarch_es_01_binance.csv")

write.csv(eth_msgarch_var_05_binance,"eth_msgarch_var_05_binance.csv")
write.csv(eth_msgarch_es_05_binance,"eth_msgarch_es_05_binance.csv")

write.csv(eth_msgarch_var_025_binance,"eth_msgarch_var_025_binance.csv")
write.csv(eth_msgarch_es_025_binance,"eth_msgarch_es_025_binance.csv")

eth_msgarch_var_01_binance= read.csv("eth_msgarch_var_01_binance.csv")
eth_msgarch_var_01_binance= eth_msgarch_var_01_binance[,-1]

eth_msgarch_vol_binance= read.csv("eth_msgarch_vol_binance.csv")
eth_msgarch_vol_binance= eth_msgarch_vol_binance[,-1]

eth_msgarch_es_01_binance= read.csv("eth_msgarch_es_01_binance.csv")
eth_msgarch_es_01_binance= eth_msgarch_es_01_binance[,-1]

eth_msgarch_var_05_binance= read.csv("eth_msgarch_var_05_binance.csv")
eth_msgarch_var_05_binance= eth_msgarch_var_05_binance[,-1]
eth_msgarch_es_05_binance= read.csv("eth_msgarch_es_05_binance.csv")
eth_msgarch_es_05_binance= eth_msgarch_es_05_binance[,-1]

eth_msgarch_var_025_binance= read.csv("eth_msgarch_var_025_binance.csv")
eth_msgarch_var_025_binance= eth_msgarch_var_025_binance[,-1]
eth_msgarch_es_025_binance= read.csv("eth_msgarch_es_025_binance.csv")
eth_msgarch_es_025_binance= eth_msgarch_es_025_binance[,-1]


# C-RV Model with normal copula,scaling factor estimation
m1= 1000
insam_est_binance_eth= c(rep(0,m1))

train_fast_sam= pk_proxy[1:1000]
ecdf= pobs(train_fast_sam)
N= length(train_fast_sam)
#ecdf= (N/(N+1))*ecd
u1= ecdf[1:(N-1)] 
u2= ecdf[2:N]
u= cbind(u1,u2)
fit.ml= fitCopula(normalCopula(), u, method="mpl")

n <- 1000
for (q in 1:1000) {
  #q = 1
  u_1 <- ecdf[q]
  U <- cCopula(cbind(u_1, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U[,2]
  my_ecdf_inv_1 = c(rep(0,n))
  for(i in 1:n){
    my_ecdf_inv_1[i] <- quantile(train_fast_sam, probs = rea[i])
  }
  insam_est_binance_eth[q]= mean(my_ecdf_inv_1)
}

insam_est_binance_eth= insam_est_binance_eth
write.csv(insam_est_binance_eth,"insam_est_binance_eth.csv")
insam_est_binance_eth= read.csv("insam_est_binance_eth.csv")
insam_est_binance_eth= insam_est_binance_eth[,-1]
scaling_factor_binance_eth= var(log_ret[1:1000])/mean(insam_est_binance_eth)


# C-RVN model forecast
nor_mean_1= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  train_fast1= pk_proxy[start:end]
  ecdf= pobs(train_fast1)
  N= length(train_fast1)
  #ecdf= (N/(N+1))*ecd
  u1= ecdf[1:(N-1)] 
  u2= ecdf[2:N]
  u= cbind(u1,u2)
  fit.ml= fitCopula(normalCopula(), u, method="mpl")
  n <- 1000
  u_1 <- ecdf[1000]
  U <- cCopula(cbind(u_1, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U[,2]
  my_ecdf_inv_1 = c(rep(0,n))
  for(i in 1:n){
    my_ecdf_inv_1[i] <- quantile(train_fast1, probs = rea[i])
  }
  nor_mean_1[q]= mean(my_ecdf_inv_1)
}

etherum_vol_nor= (nor_mean_1)
write.csv(etherum_vol_nor,"etherum_vol_nor_binance.csv")
etherum_vol_nor_binance= read.csv("etherum_vol_nor_binance.csv")
etherum_vol_nor_binance= etherum_vol_nor_binance[,-1]
adj_eth_vol_nor_binance= scaling_factor_binance_eth*(etherum_vol_nor_binance)
adj_eth_vol_nor_binance= sqrt(adj_eth_vol_nor_binance)


# C-RVt model, scaling factor for t-copula

m1= 1000
eth_insam_est_t_binance= c(rep(0,m1))

train_fast_sam= pk_proxy[1:1000]
ecdf= pobs(train_fast_sam)
N= length(train_fast_sam)
#ecdf= (N/(N+1))*ecd
u1_t= ecdf[1:(N-1)] 
u2_t= ecdf[2:N]
u_t= cbind(u1_t,u2_t)
fit.ml= fitCopula(tCopula(), u_t, method="mpl",estimate.variance=FALSE)

n <- 1000
for (q in 1:1000) {
  
  u_t <- ecdf[q]
  U_t <- cCopula(cbind(u_t, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U_t[,2]
  my_ecdf_inv_t = c(rep(0,n))
  for(i in 1:n){
    my_ecdf_inv_t[i] <- quantile(train_fast_sam, probs = rea[i])
  }
  eth_insam_est_t_binance[q]= mean(my_ecdf_inv_t)
}

eth_insam_est_t_binance= eth_insam_est_t_binance
write.csv(eth_insam_est_t_binance,"eth_insam_est_t_binance.csv")
eth_insam_est_t_binance= read.csv("eth_insam_est_t_binance.csv")
eth_insam_est_t_binance= eth_insam_est_t_binance[,-1]
eth_scaling_factor_t_binance= var(log_ret[1:1000])/mean(eth_insam_est_t_binance)

# C-RVt model forecast
t_mean_1= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  train_fast1= pk_proxy[start:end]
  ecdf= pobs(train_fast1)
  N= length(train_fast1)
  u1= ecdf[1:(N-1)]
  u2= ecdf[2:N]
  u= cbind(u1,u2)
  fit.ml= fitCopula(tCopula(), u, method="mpl",estimate.variance=FALSE)
  
  n <- 1000
  u_1 <- ecdf[1000]
  U <- cCopula(cbind(u_1, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U[,2]
  my_ecdf_inv_1 = c(rep(0,n))
  for(i in 1:n){
    my_ecdf_inv_1[i] <- quantile(train_fast1, probs = rea[i])
  }
  t_mean_1[q]= mean(my_ecdf_inv_1)
}

etherum_t_vol= (t_mean_1)
write.csv(etherum_t_vol,"etherum_t_vol_binance.csv")
etherum_t_vol= read.csv("etherum_t_vol_binance.csv")
etherum_t_vol= etherum_t_vol[,-1]
etherum_t_vol= sqrt(eth_scaling_factor_t_binance*(etherum_t_vol))

# GAS model

gas_var_01= c(rep(0,m))
gas_es_01= c(rep(0,m))
gas_var_05= c(rep(0,m))
gas_es_05= c(rep(0,m))
gas_var_025= c(rep(0,m))
gas_es_025= c(rep(0,m))
gas_vol= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  
  
  GASSpec <- UniGASSpec(Dist = "sstd", GASPar = list(scale = TRUE))
  Fit <- UniGASFit(GASSpec,res_garch)
  Forecast <- UniGASFor(Fit, H = 1)
  gas_var_for_01 <- quantile(Fit, probs = 0.01)
  gas_es_for_01 <- ES(Fit, probs = 0.01)
  gas_var_for_025 <- quantile(Fit, probs = 0.025)
  gas_es_for_025 <- ES(Fit, probs = 0.025)
  gas_var_for_05 <- quantile(Fit, probs = 0.05)
  gas_es_for_05 <- ES(Fit, probs = 0.05)
  
  gas_var_01[q]= as.numeric(tail(gas_var_for_01,1))
  gas_var_025[q]= as.numeric(tail(gas_var_for_025,1))
  gas_var_05[q]= as.numeric(tail(gas_var_for_05,1))
  gas_es_01[q]= as.numeric(tail(gas_es_for_01,1))
  gas_es_025[q]= as.numeric(tail(gas_es_for_025,1))
  gas_es_05[q]= as.numeric(tail(gas_es_for_05,1))
}

eth_gas_var_01=  gas_var_01
eth_gas_es_01=  gas_es_01
eth_gas_var_05=  gas_var_05
eth_gas_es_05=  gas_es_05
eth_gas_var_025=  gas_var_025
eth_gas_es_025=  gas_es_025

write.csv(eth_gas_var_01,"eth_gas_var_01.csv")
write.csv(eth_gas_es_01,"eth_gas_es_01.csv")
write.csv(eth_gas_var_05,"eth_gas_var_05.csv")
write.csv(eth_gas_es_05,"eth_gas_es_05.csv")
write.csv(eth_gas_var_025,"eth_gas_var_025.csv")
write.csv(eth_gas_es_025,"eth_gas_es_025.csv")


eth_gas_var_01= read.csv("eth_gas_var_01.csv")
eth_gas_var_01= eth_gas_var_01[,-1]
eth_gas_es_01= read.csv("eth_gas_es_01.csv")
eth_gas_es_01= eth_gas_es_01[,-1]

eth_gas_var_05= read.csv("eth_gas_var_05.csv")
eth_gas_var_05= eth_gas_var_05[,-1]
eth_gas_es_05= read.csv("eth_gas_es_05.csv")
eth_gas_es_05= eth_gas_es_05[,-1]

eth_gas_var_025= read.csv("eth_gas_var_025.csv")
eth_gas_var_025= eth_gas_var_025[,-1]
eth_gas_es_025= read.csv("eth_gas_es_025.csv")
eth_gas_es_025= eth_gas_es_025[,-1]

# RooBGARCH Model

library(RobGARCHBoot)

VaR_rob_01= c(rep(0,m))
VaR_rob_05= c(rep(0,m))
VaR_rob_025= c(rep(0,m))

ES_rob_01= c(rep(0,m))
ES_rob_05= c(rep(0,m))
ES_rob_025= c(rep(0,m))

for(q in 1:m){
  
  start=q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  res_garch= fit$residuals
  
  boot = RobGARCHBoot(res_garch, n.boot = 1000, n.ahead = 1)
  return= boot[[1]]
  VaR_rob_01[q] = quantile(boot[[1]], prob = 0.01)
  VaR_rob_05[q] = quantile(boot[[1]], prob = 0.05)
  VaR_rob_025[q] = quantile(boot[[1]], prob = 0.025)
  
  ES_rob_01[q]= mean(return[return<VaR_rob_01[q]])
  ES_rob_05[q]= mean(return[return<VaR_rob_05[q]])
  ES_rob_025[q]= mean(return[return<VaR_rob_025[q]])
  
  
  
}

eth_rob_garch_var_01= (VaR_rob_01)
eth_rob_garch_var_05= (VaR_rob_05)
eth_rob_garch_var_025= (VaR_rob_025)

eth_rob_garch_ES_01= (ES_rob_01)
eth_rob_garch_ES_05= (ES_rob_05)
eth_rob_garch_ES_025= (ES_rob_025)

write.csv(eth_rob_garch_var_01,"eth_rob_garch_var_01.csv")
write.csv(eth_rob_garch_var_05,"eth_rob_garch_var_05.csv")
write.csv(eth_rob_garch_var_025,"eth_rob_garch_var_025.csv")

write.csv(eth_rob_garch_ES_01,"eth_rob_garch_ES_01.csv")
write.csv(eth_rob_garch_ES_05,"eth_rob_garch_ES_05.csv")
write.csv(eth_rob_garch_ES_025,"eth_rob_garch_ES_025.csv")

eth_sq_vol_rob= sq_vol_rob
write.csv(eth_sq_vol_rob,"eth_sq_vol_rob.csv")

eth_rob_garch_var_01= read.csv("eth_rob_garch_var_01.csv")
eth_rob_garch_var_01= eth_rob_garch_var_01[,-1]
eth_rob_garch_var_025= read.csv("eth_rob_garch_var_025.csv")
eth_rob_garch_var_025= eth_rob_garch_var_025[,-1]
eth_rob_garch_var_05= read.csv("eth_rob_garch_var_05.csv")
eth_rob_garch_var_05= eth_rob_garch_var_05[,-1]

eth_rob_garch_es_01= read.csv("eth_rob_garch_es_01.csv")
eth_rob_garch_es_01= eth_rob_garch_es_01[,-1]
eth_rob_garch_es_025= read.csv("eth_rob_garch_es_025.csv")
eth_rob_garch_es_025= eth_rob_garch_es_025[,-1]
eth_rob_garch_es_05= read.csv("eth_rob_garch_es_05.csv")
eth_rob_garch_es_05= eth_rob_garch_es_05[,-1]


# CARR Model in scaling factor estimation

in_sam_ran= range[1:1000]
spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
e_garch_range <- ugarchfit(in_sam_ran, spec=spec.gjrGARCH_range)
est_range= (e_garch_range@fit$sigma)^2
scaling_carr_binance= sd(log_ret[1:1000])/mean(est_range)

# CARR Model forecast
carr_vol_binance= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  train_fast1= range[start:end]
  spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
  e_garch_range <- ugarchfit(train_fast1, spec=spec.gjrGARCH_range)
  predi_garch_range= ugarchforecast(e_garch_range,n.ahead=1,data=train_fast1)
  fore= (predi_garch_range@forecast)
  carr_vol_binance[q]= as.numeric(fore$sigma)
}

etherum_carr_vol_binance= (carr_vol_binance)
write.csv(etherum_carr_vol_binance,"etherum_carr_vol_binance.csv")
etherum_carr_vol_binance= read.csv("etherum_carr_vol_binance.csv")
etherum_carr_vol_binance= etherum_carr_vol_binance[,-1]
etherum_carr_vol_binance= (etherum_carr_vol_binance)^2
adj_etherum_carr_vol_binance= scaling_carr_binance*(etherum_carr_vol_binance)

# GARCH-ged model

pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

# GARCH model forecast with sstd distribution

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  fit_mean_garch[q]= as.numeric(pre_mean_garch$pred)
  
  spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH)
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_garch_fast1[q]= as.numeric(fore$sigmaFor)
  
  
}

etherum_garch_mean_ged= fit_mean_garch
etherum_garch_ged= pre_garch_fast1

write.csv(etherum_garch_mean_ged,"etherum_garch_mean_ged.csv")
etherum_garch_mean= read.csv("etherum_garch_mean_ged.csv")
etherum_garch_mean= etherum_garch_mean[,-1]
write.csv(etherum_garch_ged,"etherum_garch_ged.csv")
etherum_garch_vol_ged= read.csv("etherum_garch_ged.csv")
etherum_garch_vol_ged= etherum_garch_vol_ged[,-1]

# NAGARCH-ged Model

pre_nagarch_binance= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  
  spec.gjrGARCH = ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(1,1),submodel= "NAGARCH"), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH,solver = "hybrid")
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_nagarch_binance[q]= as.numeric(fore$sigmaFor)
  
}

bid_nagarch_binance_eth_ged= pre_nagarch_binance
write.csv(bid_nagarch_binance_eth_ged,"bid_nagarch_binance_eth_ged.csv")
nagarch_binance_eth_ged= read.csv("bid_nagarch_binance_eth_ged.csv")
nagarch_binance_eth_ged= nagarch_binance_eth_ged[,-1]

#GJR-GARCH-ged Model

pre_gjrgarch_fast1=c(rep(0,m))
for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  res_gjrgarch= fit$residuals
  spec.garch = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_gjrgarch <- ugarchfit(res_gjrgarch, spec=spec.garch)
  e_forecast= ugarchforecast(e_gjrgarch,n.ahead=1,data=res_gjrgarch)
  fore= e_forecast@forecast
  pre_gjrgarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

eth_gjrgarch_ged= pre_gjrgarch_fast1
write.csv(eth_gjrgarch_ged,"eth_gjrgarch_ged.csv")
eth_gjrgarch_ged= read.csv("eth_gjrgarch_ged.csv")
eth_gjrgarch_ged= eth_gjrgarch_ged[,-1]

# VaR and ES estimation for different Models and backtesting


# VaR with sstd

fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 3,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="sstd")
e_garch_sstd<- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)



VaR_sstd_garch =  (etherum_garch_vol)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                            shape=coef(e_garch_sstd)["shape"])

VaR_sstd_nor_cop =  (adj_eth_vol_nor_binance)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                                    shape=coef(e_garch_sstd)["shape"])


VaR_sstd_t_cop =  (etherum_t_vol)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                        shape=coef(e_garch_sstd)["shape"])



VaR_sstd_carr =   (adj_etherum_carr_vol_binance)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                                       shape=coef(e_garch_sstd)["shape"])

VaR_sstd_nagarch =  (nagarch_binance_eth)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                                shape=coef(e_garch_sstd)["shape"])
VaR_sstd_gjrgarch =  (eth_gjrgarch)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],shape=coef(e_garch_sstd)["shape"])



# VaR with ged

fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 3,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="ged")
e_garch <- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)


VaR_ged_garch =  (etherum_garch_vol)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                           shape=coef(e_garch)["shape"])

VaR_ged_nor_cop =  (adj_eth_vol_nor_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                                   shape=coef(e_garch)["shape"])


VaR_ged_t_cop =  (etherum_t_vol)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                       shape=coef(e_garch)["shape"])



VaR_ged_carr =   (adj_etherum_carr_vol_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                                      shape=coef(e_garch)["shape"])

VaR_ged_nagarch =  (nagarch_binance_eth)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                               
                                               shape=coef(e_garch)["shape"])
VaR_ged_gjrgarch =  (eth_gjrgarch)*qdist("ged", p=0.05, mu = 0, sigma = 1,shape=coef(e_garch)["shape"])

test= log_ret[1001:2343]
actual= test-etherum_garch_mean

backtest_sstd_nor_cop= BacktestVaR(data = actual, VaR = VaR_sstd_nor_cop, alpha = 0.05)

backtest_sstd_garch= BacktestVaR(data = actual, VaR = VaR_sstd_garch, alpha = 0.05)

backtest_sstd_carr= BacktestVaR(data = actual, VaR = VaR_sstd_carr, alpha = 0.05)

backtest_sstd_t_cop= BacktestVaR(data = actual, VaR = VaR_sstd_t_cop, alpha = 0.05)

backtest_robgarch= BacktestVaR(data = actual, VaR = eth_rob_garch_var_05, alpha = 0.05)

backtest_sstd_nagarch= BacktestVaR(data = actual, VaR = VaR_sstd_nagarch, alpha = 0.05)

backtest_ged_nor_cop= BacktestVaR(data = actual, VaR = VaR_ged_nor_cop, alpha = 0.05)

backtest_ged_garch= BacktestVaR(data = actual, VaR = VaR_ged_garch, alpha = 0.05)

backtest_ged_carr= BacktestVaR(data = actual, VaR = VaR_ged_carr, alpha = 0.05)

backtest_ged_t_cop= BacktestVaR(data = actual, VaR = VaR_ged_t_cop, alpha = 0.05)
backtest_ged_nagarch= BacktestVaR(data = actual, VaR = VaR_ged_nagarch, alpha = 0.05)
backtest_ged_gjrgarch= BacktestVaR(data = actual, VaR = VaR_ged_gjrgarch, alpha = 0.05)

backtest_msgarch= BacktestVaR(data = actual, VaR = eth_msgarch_var_05_binance, alpha = 0.05)

backtest_gas= BacktestVaR(data = actual, VaR = eth_gas_var_05, alpha = 0.05)

backtest_sstd_gjrgarch= BacktestVaR(data = actual, VaR = VaR_sstd_gjrgarch, alpha = 0.05)


# UC Test
uc_ged_nor_cop= (backtest_ged_nor_cop$LRuc["Pvalue" ])
uc_ged_t_cop= (backtest_ged_t_cop$LRuc["Pvalue" ])
uc_ged_garch= backtest_ged_garch$LRuc["Pvalue" ]
uc_ged_nagarch= backtest_ged_nagarch$LRuc["Pvalue" ]
uc_ged_gjrgarch= backtest_ged_gjrgarch$LRuc["Pvalue" ]
uc_ged_carr= backtest_ged_carr$LRuc["Pvalue" ]
uc_msgarch= backtest_msgarch$LRuc["Pvalue" ]
uc_robgarch= backtest_robgarch$LRuc["Pvalue" ]
uc_gas= backtest_gas$LRuc["Pvalue" ]
uc_sstd_nor_cop= (backtest_sstd_nor_cop$LRuc["Pvalue" ])
uc_sstd_t_cop= (backtest_sstd_t_cop$LRuc["Pvalue" ])
uc_sstd_garch= backtest_sstd_garch$LRuc["Pvalue" ]
uc_sstd_nagarch= backtest_sstd_nagarch$LRuc["Pvalue" ]
uc_sstd_gjrgarch= backtest_sstd_gjrgarch$LRuc["Pvalue" ]
uc_sstd_carr= backtest_sstd_carr$LRuc["Pvalue" ]

# CC Test
cc_ged_nor_cop= (backtest_ged_nor_cop$LRcc["Pvalue" ])
cc_ged_t_cop= (backtest_ged_t_cop$LRcc["Pvalue" ])
cc_ged_garch= backtest_ged_garch$LRcc["Pvalue" ]
cc_ged_nagarch= backtest_ged_nagarch$LRcc["Pvalue" ]
cc_ged_gjrgarch= backtest_ged_gjrgarch$LRcc["Pvalue" ]
cc_ged_carr= backtest_ged_carr$LRcc["Pvalue" ]
cc_msgarch= backtest_msgarch$LRcc["Pvalue" ]
cc_robgarch= backtest_robgarch$LRcc["Pvalue" ]
cc_gas= backtest_gas$LRcc["Pvalue" ]
cc_sstd_nor_cop= (backtest_sstd_nor_cop$LRcc["Pvalue" ])
cc_sstd_t_cop= (backtest_sstd_t_cop$LRcc["Pvalue" ])
cc_sstd_garch= backtest_sstd_garch$LRcc["Pvalue" ]
cc_sstd_nagarch= backtest_sstd_nagarch$LRcc["Pvalue" ]
cc_sstd_gjrgarch= backtest_sstd_gjrgarch$LRcc["Pvalue" ]
cc_sstd_carr= backtest_sstd_carr$LRcc["Pvalue" ]

# DQ Test
DQ_ged_nor_cop= as.numeric(backtest_ged_nor_cop$DQ$pvalue)
DQ_ged_t_cop= as.numeric(backtest_ged_t_cop$DQ$pvalue)
DQ_ged_garch=as.numeric( backtest_ged_garch$DQ$pvalue)
DQ_ged_nagarch=as.numeric( backtest_ged_nagarch$DQ$pvalue)
DQ_ged_gjrgarch=as.numeric( backtest_ged_gjrgarch$DQ$pvalue)
DQ_ged_carr=as.numeric( backtest_ged_carr$DQ$pvalue)
DQ_msgarch=as.numeric( backtest_msgarch$DQ$pvalue)
DQ_robgarch=as.numeric( backtest_robgarch$DQ$pvalue)
DQ_gas=as.numeric( backtest_gas$DQ$pvalue)
DQ_sstd_nor_cop= as.numeric(backtest_sstd_nor_cop$DQ$pvalue)
DQ_sstd_t_cop= as.numeric(backtest_sstd_t_cop$DQ$pvalue)
DQ_sstd_garch=as.numeric( backtest_sstd_garch$DQ$pvalue)
DQ_sstd_nagarch=as.numeric( backtest_sstd_nagarch$DQ$pvalue)
DQ_sstd_gjrgarch=as.numeric( backtest_sstd_gjrgarch$DQ$pvalue)
DQ_sstd_carr=as.numeric( backtest_sstd_carr$DQ$pvalue)

# QL Value
QL_ged_nor_cop= as.numeric(backtest_ged_nor_cop$Loss$Loss)
QL_ged_t_cop= as.numeric(backtest_ged_t_cop$Loss$Loss)
QL_ged_garch=as.numeric( backtest_ged_garch$Loss$Loss)
QL_ged_nagarch=as.numeric( backtest_ged_nagarch$Loss$Loss)
QL_ged_gjrgarch=as.numeric( backtest_ged_gjrgarch$Loss$Loss)
QL_ged_carr=as.numeric( backtest_ged_carr$Loss$Loss)
QL_msgarch=as.numeric( backtest_msgarch$Loss$Loss)
QL_robgarch=as.numeric( backtest_robgarch$Loss$Loss)
QL_gas=as.numeric( backtest_gas$Loss$Loss)
QL_sstd_nor_cop= as.numeric(backtest_sstd_nor_cop$Loss$Loss)
QL_sstd_t_cop= as.numeric(backtest_sstd_t_cop$Loss$Loss)
QL_sstd_garch=as.numeric( backtest_sstd_garch$Loss$Loss)
QL_sstd_nagarch=as.numeric( backtest_sstd_nagarch$Loss$Loss)
QL_sstd_gjrgarch=as.numeric( backtest_sstd_gjrgarch$Loss$Loss)
QL_sstd_carr=as.numeric( backtest_sstd_carr$Loss$Loss)

# Expected shortfall and McF Test

f1 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])

ES_sstd_garch = (etherum_garch_vol)*integrate(f1, 0, 0.05)$value/0.05

# McF test
mcf_sstd_garch= (ESTest(0.05, actual, ES_sstd_garch, VaR_sstd_garch, boot = TRUE))



f2 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_nor_vol = (adj_eth_vol_nor_binance)*integrate(f2, 0, 0.05)$value/0.05
# McF test
mcf_sstd_nor_cop= (ESTest(0.05, actual, ES_sstd_nor_vol, VaR_sstd_nor_cop, boot = TRUE))



f4 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_t_vol = (etherum_t_vol)*integrate(f4, 0, 0.05)$value/0.05
# McF test

mcf_sstd_t_cop=(ESTest(0.05, actual, ES_sstd_t_vol,VaR_sstd_t_cop, boot = TRUE))

# McF test

mcf_robgarch =(ESTest(0.05, actual, eth_rob_garch_es_05, eth_rob_garch_var_05, boot = TRUE))


f3 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_carr = (adj_etherum_carr_vol_binance)*integrate(f3, 0, 0.05)$value/0.05
# McF test
mcf_sstd_carr= (ESTest(0.05, actual, ES_sstd_carr, VaR_sstd_carr, boot = TRUE))

f5 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_nagarch = (nagarch_binance_eth)*integrate(f5, 0, 0.05)$value/0.05
# McF test
mcf_sstd_nagarch= (ESTest(0.05, actual, ES_sstd_nagarch, VaR_sstd_nagarch, boot = TRUE))

# McF test
mcf_msgarch=(ESTest(0.05, actual, eth_msgarch_es_05_binance, eth_msgarch_var_05_binance, boot = TRUE))

mcf_gas =(ESTest(0.05, test, eth_gas_es_05, eth_gas_var_05, boot = TRUE))

f6 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])

ES_sstd_gjrgarch = (eth_gjrgarch)*integrate(f6, 0, 0.05)$value/0.05

# McF test
mcf_sstd_gjrgarch=(ESTest(0.05, actual, ES_sstd_gjrgarch, VaR_sstd_gjrgarch, boot = TRUE))

f1_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])

ES_ged_garch = (etherum_garch_vol_ged)*integrate(f1_ged, 0, 0.05)$value/0.05

# McF test
mcf_ged_garch= (ESTest(0.05, actual, ES_ged_garch, VaR_ged_garch, boot = TRUE))



f2_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_nor_vol = (adj_eth_vol_nor_binance)*integrate(f2_ged, 0, 0.05)$value/0.05
# McF test
mcf_ged_nor_cop= (ESTest(0.05, actual, ES_ged_nor_vol, VaR_ged_nor_cop, boot = TRUE))



f4_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_t_vol = (etherum_t_vol)*integrate(f4_ged, 0, 0.05)$value/0.05
# McF test

mcf_ged_t_cop=(ESTest(0.05, actual, ES_ged_t_vol,VaR_ged_t_cop, boot = TRUE))

f3_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_carr = (adj_etherum_carr_vol_binance)*integrate(f3_ged, 0, 0.05)$value/0.05
# McF test
mcf_ged_carr= (ESTest(0.05, actual, ES_ged_carr, VaR_ged_carr, boot = TRUE))

f5_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_nagarch = (nagarch_binance_eth_ged)*integrate(f5_ged, 0, 0.05)$value/0.05
# McF test
mcf_ged_nagarch= (ESTest(0.05, actual, ES_ged_nagarch, VaR_ged_nagarch, boot = TRUE))

f6_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])

ES_ged_gjrgarch = (eth_gjrgarch_ged)*integrate(f6_ged, 0, 0.05)$value/0.05

# McF test
mcf_ged_gjrgarch=(ESTest(0.05, actual, ES_ged_gjrgarch, VaR_ged_gjrgarch, boot = TRUE))


# Percentage of Hits

hits_nor_ged_cop= (mcf_ged_nor_cop$actual.exceed*100)/1343
hits_t_ged_cop= (mcf_ged_t_cop$actual.exceed*100)/1343
hits_ged_garch= (mcf_ged_garch$actual.exceed*100)/1343
hits_ged_nagarch= (mcf_ged_nagarch$actual.exceed*100)/1343
hits_ged_gjrgarch= (mcf_ged_gjrgarch$actual.exceed*100)/1343
hits_ged_carr= (mcf_ged_carr$actual.exceed*100)/1343
hits_msgarch= (mcf_msgarch$actual.exceed*100)/1343
hits_robgarch= (mcf_robgarch$actual.exceed*100)/1343
hits_gas= (mcf_gas$actual.exceed*100)/1343
hits_nor_sstd_cop= (mcf_sstd_nor_cop$actual.exceed*100)/1343
hits_t_sstd_cop= (mcf_sstd_t_cop$actual.exceed*100)/1343
hits_sstd_garch= (mcf_sstd_garch$actual.exceed*100)/1343
hits_sstd_nagarch= (mcf_sstd_nagarch$actual.exceed*100)/1343
hits_sstd_gjrgarch= (mcf_sstd_gjrgarch$actual.exceed*100)/1343
hits_sstd_carr= (mcf_sstd_carr$actual.exceed*100)/1343


# McF Test

McF_p_ged_nor_cop= mcf_ged_nor_cop$p.value
McF_p_ged_t_cop= mcf_ged_t_cop$p.value
McF_p_ged_garch= mcf_ged_garch$p.value
McF_p_ged_nagarch= mcf_ged_nagarch$p.value
McF_p_ged_gjrgarch= mcf_ged_gjrgarch$p.value
McF_p_ged_carr= mcf_ged_carr$p.value
McF_p_msgarch= mcf_msgarch$p.value
McF_p_robgarch= mcf_robgarch$p.value
McF_p_gas= mcf_gas$p.value
McF_p_sstd_nor_cop= mcf_sstd_nor_cop$p.value
McF_p_sstd_t_cop= mcf_sstd_t_cop$p.value
McF_p_sstd_garch= mcf_sstd_garch$p.value
McF_p_sstd_nagarch= mcf_sstd_nagarch$p.value
McF_p_sstd_gjrgarch= mcf_sstd_gjrgarch$p.value
McF_p_sstd_carr= mcf_sstd_carr$p.value

# NF Test
NF_sstd_nor_cop= cc_backtest(actual, VaR_sstd_nor_cop, ES_sstd_nor_vol, alpha=0.05)
NF_sstd_t_cop= cc_backtest(actual, VaR_sstd_t_cop, ES_sstd_t_vol, alpha=0.05)
NF_sstd_garch= cc_backtest(actual, VaR_sstd_garch, ES_sstd_garch, alpha = 0.05)
NF_sstd_gjrgarch= cc_backtest(actual, VaR_sstd_gjrgarch, ES_sstd_gjrgarch,alpha = 0.05)
NF_sstd_nagarch= cc_backtest(actual, VaR_sstd_nagarch, ES_sstd_nagarch,alpha = 0.05)
NF_sstd_carr= cc_backtest(actual, VaR_sstd_carr, ES_sstd_carr, alpha=0.05)

NF_msgarch= cc_backtest(actual, eth_msgarch_var_05_binance, eth_msgarch_es_05_binance, alpha=0.05)
NF_gas= cc_backtest(actual, eth_gas_var_05, eth_gas_es_05, alpha=0.05)
NF_robgarch= cc_backtest(actual, eth_rob_garch_var_05, eth_rob_garch_es_05, alpha=0.05)
NF_ged_nor_cop= cc_backtest(actual, VaR_ged_nor_cop, ES_ged_nor_vol, alpha=0.05)
NF_ged_t_cop= cc_backtest(actual, VaR_ged_t_cop, ES_ged_t_vol, alpha=0.05)
NF_ged_garch= cc_backtest(actual, VaR_ged_garch, ES_ged_garch, alpha = 0.05)
NF_ged_gjrgarch= cc_backtest(actual, VaR_ged_gjrgarch, ES_ged_gjrgarch,alpha = 0.05)
NF_ged_nagarch= cc_backtest(actual, VaR_ged_nagarch, ES_ged_nagarch,alpha = 0.05)
NF_ged_carr= cc_backtest(actual, VaR_ged_carr, ES_ged_carr, alpha=0.05)

#NF p-value

NF_p_ged_nor_cop= NF_ged_nor_cop$pvalue_twosided_simple
NF_p_ged_t_cop= NF_ged_t_cop$pvalue_twosided_simple
NF_p_ged_garch= NF_ged_garch$pvalue_twosided_simple
NF_p_ged_nagarch= NF_ged_nagarch$pvalue_twosided_simple
NF_p_ged_gjrgarch= NF_ged_gjrgarch$pvalue_twosided_simple
NF_p_ged_carr= NF_ged_carr$pvalue_twosided_simple
NF_p_msgarch= NF_msgarch$pvalue_twosided_simple
NF_p_robgarch= NF_robgarch$pvalue_twosided_simple
NF_p_gas= NF_gas$pvalue_twosided_simple
NF_p_sstd_nor_cop= NF_sstd_nor_cop$pvalue_twosided_simple
NF_p_sstd_t_cop= NF_sstd_t_cop$pvalue_twosided_simple
NF_p_sstd_garch= NF_sstd_garch$pvalue_twosided_simple
NF_p_sstd_nagarch= NF_sstd_nagarch$pvalue_twosided_simple
NF_p_sstd_gjrgarch= NF_sstd_gjrgarch$pvalue_twosided_simple
NF_p_sstd_carr= NF_sstd_carr$pvalue_twosided_simple

# BD Test

BD_sstd_nor_cop= esr_backtest(actual, VaR_sstd_nor_cop, ES_sstd_nor_vol, alpha=0.05,version = 1)
BD_sstd_t_cop= esr_backtest(actual, VaR_sstd_t_cop, ES_sstd_t_vol, alpha=0.05,version = 1)
BD_sstd_garch= esr_backtest(actual, VaR_sstd_garch, ES_sstd_garch, alpha = 0.05,version = 1)
BD_sstd_gjrgarch= esr_backtest(actual, VaR_sstd_gjrgarch, ES_sstd_gjrgarch,alpha = 0.05,version = 1)
BD_sstd_nagarch= esr_backtest(actual, VaR_sstd_nagarch, ES_sstd_nagarch,alpha = 0.05,version = 1)
BD_sstd_carr= esr_backtest(actual, VaR_sstd_carr, ES_sstd_carr, alpha=0.05,version = 1)
BD_msgarch= esr_backtest(actual, eth_msgarch_var_05_binance, eth_msgarch_es_05_binance, alpha=0.05,version = 1)
BD_gas= esr_backtest(actual, eth_gas_var_05, eth_gas_es_05, alpha=0.05,version = 1)
BD_robgarch= esr_backtest(actual, eth_rob_garch_var_05, eth_rob_garch_es_05, alpha=0.05,version = 1)
BD_ged_nor_cop= esr_backtest(actual, VaR_ged_nor_cop, ES_ged_nor_vol, alpha=0.05,version = 1)
BD_ged_t_cop= esr_backtest(actual, VaR_ged_t_cop, ES_ged_t_vol, alpha=0.05,version = 1)
BD_ged_garch= esr_backtest(actual, VaR_ged_garch, ES_ged_garch, alpha = 0.05,version = 1)
BD_ged_gjrgarch= esr_backtest(actual, VaR_ged_gjrgarch, ES_ged_gjrgarch,alpha = 0.05,version = 1)
BD_ged_nagarch= esr_backtest(actual, VaR_ged_nagarch, ES_ged_nagarch,alpha = 0.05,version = 1)
BD_ged_carr= esr_backtest(actual, VaR_ged_carr, ES_ged_carr, alpha=0.05,version = 1)

# BD Test p-value

BD_p_ged_nor_cop= BD_ged_nor_cop$pvalue_twosided_asymptotic
BD_p_ged_t_cop= BD_ged_t_cop$pvalue_twosided_asymptotic
BD_p_ged_garch= BD_ged_garch$pvalue_twosided_asymptotic
BD_p_ged_nagarch= BD_ged_nagarch$pvalue_twosided_asymptotic
BD_p_ged_gjrgarch= BD_ged_gjrgarch$pvalue_twosided_asymptotic
BD_p_ged_carr= BD_ged_carr$pvalue_twosided_asymptotic
BD_p_msgarch= BD_msgarch$pvalue_twosided_asymptotic
BD_p_robgarch= BD_robgarch$pvalue_twosided_asymptotic
BD_p_gas= BD_gas$pvalue_twosided_asymptotic
BD_p_sstd_nor_cop= BD_sstd_nor_cop$pvalue_twosided_asymptotic
BD_p_sstd_t_cop= BD_sstd_t_cop$pvalue_twosided_asymptotic
BD_p_sstd_garch= BD_sstd_garch$pvalue_twosided_asymptotic
BD_p_sstd_nagarch= BD_sstd_nagarch$pvalue_twosided_asymptotic
BD_p_sstd_gjrgarch= BD_sstd_gjrgarch$pvalue_twosided_asymptotic
BD_p_sstd_carr= BD_sstd_carr$pvalue_twosided_asymptotic

# FZL Function
FZL_sstd_nor_cop= mean(FZLoss(actual, VaR_sstd_nor_cop, ES_sstd_nor_vol, alpha=0.05))
FZL_sstd_t_cop= mean(FZLoss(actual, VaR_sstd_t_cop, ES_sstd_t_vol, alpha=0.05))
FZL_sstd_garch= mean(FZLoss(actual, VaR_sstd_garch, ES_sstd_garch, alpha = 0.05))
FZL_sstd_gjrgarch= mean(FZLoss(actual, VaR_sstd_gjrgarch, ES_sstd_gjrgarch,alpha = 0.05))
FZL_sstd_nagarch= mean(FZLoss(actual, VaR_sstd_nagarch, ES_sstd_nagarch,alpha = 0.05))
FZL_sstd_carr= mean(FZLoss(actual, VaR_sstd_carr, ES_sstd_carr, alpha=0.05))
FZL_msgarch= mean(FZLoss(actual, eth_msgarch_var_05_binance, eth_msgarch_es_05_binance, alpha=0.05))
FZL_gas= mean(FZLoss(actual, eth_gas_var_05, eth_gas_es_05, alpha=0.05))
FZL_robgarch= mean(FZLoss(actual, eth_rob_garch_var_05, eth_rob_garch_es_05, alpha=0.05))
FZL_ged_nor_cop= mean(FZLoss(actual, VaR_ged_nor_cop, ES_ged_nor_vol, alpha=0.05))
FZL_ged_t_cop= mean(FZLoss(actual, VaR_ged_t_cop, ES_ged_t_vol, alpha=0.05))
FZL_ged_garch= mean(FZLoss(actual, VaR_ged_garch, ES_ged_garch, alpha = 0.05))
FZL_ged_gjrgarch= mean(FZLoss(actual, VaR_ged_gjrgarch, ES_ged_gjrgarch,alpha = 0.05))
FZL_ged_nagarch= mean(FZLoss(actual, VaR_ged_nagarch, ES_ged_nagarch,alpha = 0.05))
FZL_ged_carr= mean(FZLoss(actual, VaR_ged_carr, ES_ged_carr, alpha=0.05))

Hits= rbind(hits_nor_ged_cop,hits_t_ged_cop,hits_ged_garch,hits_ged_nagarch,hits_ged_gjrgarch,hits_ged_carr,hits_msgarch,hits_robgarch,hits_gas,hits_nor_sstd_cop,hits_t_sstd_cop,hits_sstd_garch,hits_sstd_nagarch,hits_sstd_gjrgarch,hits_sstd_carr)
UC= rbind(uc_ged_nor_cop,uc_ged_t_cop,uc_ged_garch,uc_ged_nagarch,uc_ged_gjrgarch,uc_ged_carr,uc_msgarch,uc_robgarch,uc_gas,uc_sstd_nor_cop,uc_sstd_t_cop,uc_sstd_garch,uc_sstd_nagarch,uc_sstd_gjrgarch,uc_sstd_carr)
CC= rbind(cc_ged_nor_cop,cc_ged_t_cop,cc_ged_garch,cc_ged_nagarch,cc_ged_gjrgarch,cc_ged_carr,cc_msgarch,cc_robgarch,cc_gas,cc_sstd_nor_cop,cc_sstd_t_cop,cc_sstd_garch,cc_sstd_nagarch,cc_sstd_gjrgarch,cc_sstd_carr)
DQ= rbind(DQ_ged_nor_cop,DQ_ged_t_cop,DQ_ged_garch,DQ_ged_nagarch,DQ_ged_gjrgarch,DQ_ged_carr,DQ_msgarch,DQ_robgarch,DQ_gas,DQ_sstd_nor_cop,DQ_sstd_t_cop,DQ_sstd_garch,DQ_sstd_nagarch,DQ_sstd_gjrgarch,DQ_sstd_carr)
McF= rbind(McF_p_ged_nor_cop,McF_p_ged_t_cop,McF_p_ged_garch,McF_p_ged_nagarch,McF_p_ged_gjrgarch,McF_p_ged_carr,McF_p_msgarch,McF_p_robgarch,McF_p_gas,McF_p_sstd_nor_cop,McF_p_sstd_t_cop,McF_p_sstd_garch,McF_p_sstd_nagarch,McF_p_sstd_gjrgarch,McF_p_sstd_carr)
NF= rbind(NF_p_ged_nor_cop,NF_p_ged_t_cop,NF_p_ged_garch,NF_p_ged_nagarch,NF_p_ged_gjrgarch,NF_p_ged_carr,NF_p_msgarch,NF_p_robgarch,NF_p_gas,NF_p_sstd_nor_cop,NF_p_sstd_t_cop,NF_p_sstd_garch,NF_p_sstd_nagarch,NF_p_sstd_gjrgarch,NF_p_sstd_carr)
BD= rbind(BD_p_ged_nor_cop,BD_p_ged_t_cop,BD_p_ged_garch,BD_p_ged_nagarch,BD_p_ged_gjrgarch,BD_p_ged_carr,BD_p_msgarch,BD_p_robgarch,BD_p_gas,BD_p_sstd_nor_cop,BD_p_sstd_t_cop,BD_p_sstd_garch,BD_p_sstd_nagarch,BD_p_sstd_gjrgarch,BD_p_sstd_carr)
QL= rbind(QL_ged_nor_cop,QL_ged_t_cop,QL_ged_garch,QL_ged_nagarch,QL_ged_gjrgarch,QL_ged_carr,QL_msgarch,QL_robgarch,QL_gas,QL_sstd_nor_cop,QL_sstd_t_cop,QL_sstd_garch,QL_sstd_nagarch,QL_sstd_gjrgarch,QL_sstd_carr)
FZL= rbind(FZL_ged_nor_cop,FZL_ged_t_cop,FZL_ged_garch,FZL_ged_nagarch,FZL_ged_gjrgarch,FZL_ged_carr,FZL_msgarch,FZL_robgarch,FZL_gas,FZL_sstd_nor_cop,FZL_sstd_t_cop,FZL_sstd_garch,FZL_sstd_nagarch,FZL_sstd_gjrgarch,FZL_sstd_carr)

Table_eth= cbind(Hits,UC,CC,DQ,McF,NF,BD,QL,FZL)

# Similarly estime VaR and ES at 1% and 2.5% risk levels