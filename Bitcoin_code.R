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

data <- read.csv("Bitcoin 2017 data.csv")

l= length(data$Close)
log_ret= c(rep(0,l))

for(i in 2:l){
  log_ret[i]= log(data$Close[i]/data$Close[i-1])*100
}
log_ret= log_ret[-1]

high= data$High
low= data$Low

# Parkinson estimator
pk_proxy= {(log(high)-log(low))^2}/(4*log(2))
pk_proxy= pk_proxy[-1]

# Range estimator for CARR model
log_range= c(rep(0,l))

for(i in 1:l){
  log_range[i]= log(high[i])-log(low[(i)])
}

log_range=log_range[-1]

range= sqrt(log_range)

# Rolling window forecast for GARCH-sstd Model

m= 1343
pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

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

garch_mean= fit_mean_garch
bid_garch= pre_garch_fast1
write.csv(bid_garch,"bid_garch.csv")
write.csv(garch_mean,"garch_mean.csv")
bid_garch= read.csv("bid_garch.csv")
bid_garch= bid_garch[,-1]
bid_mean= read.csv("garch_mean.csv")
bid_mean= bid_mean[,-1]

# Rolling window forecast for NAGARCH-sstd Model


pre_nagarch_fast1= c(rep(0,m))

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
  pre_nagarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

bid_nagarch= pre_nagarch_fast1
write.csv(bid_nagarch,"bid_nagarch.csv")
bid_nagarch= read.csv("bid_nagarch.csv")
bid_nagarch= bid_nagarch[,-1]

#GJR-GARCH-sstd Model

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

bid_gjrgarch= pre_gjrgarch_fast1
write.csv(bid_gjrgarch,"bid_gjrgarch.csv")
bid_gjrgarch= read.csv("bid_gjrgarch.csv")
bid_gjrgarch= bid_gjrgarch[,-1]


#  Rolling window forecast for MSGARCH Model
msgarch_var_01= c(rep(0,m))
msgarch_es_01= c(rep(0,m))
msgarch_var_05= c(rep(0,m))
msgarch_es_05= c(rep(0,m))
msgarch_var_025= c(rep(0,m))
msgarch_es_025= c(rep(0,m))
msgarch_vol= c(rep(0,m))

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
  msgarch_vol[q]= ms_garch$vol
  
  e_forecast_01= Risk(fit.ml, alpha = 0.01, nahead = 1)
  e_forecast_05= Risk(fit.ml, alpha = 0.05, nahead = 1)
  e_forecast_025= Risk(fit.ml, alpha = 0.025, nahead = 1)
  
  msgarch_var_01[q]= e_forecast_01$VaR
  msgarch_es_01[q]= e_forecast_01$ES
  msgarch_var_05[q]= e_forecast_05$VaR
  msgarch_es_05[q]= e_forecast_05$ES
  
  msgarch_var_025[q]= e_forecast_025$VaR
  msgarch_es_025[q]= e_forecast_025$ES
  
}

bid_msgarch_var_01=  msgarch_var_01
bid_msgarch_es_01=  msgarch_es_01
bid_msgarch_var_05=  msgarch_var_05
bid_msgarch_es_05=  msgarch_es_05

bid_msgarch_var_025=  msgarch_var_025
bid_msgarch_es_025=  msgarch_es_025

msgarch_vol= msgarch_vol

write.csv(msgarch_vol,"msgarch_vol.csv")

write.csv(bid_msgarch_var_01,"bid_msgarch_var_01.csv")
write.csv(bid_msgarch_es_01,"bid_msgarch_es_01.csv")

write.csv(bid_msgarch_var_05,"bid_msgarch_var_05.csv")
write.csv(bid_msgarch_es_05,"bid_msgarch_es_05.csv")

write.csv(bid_msgarch_var_025,"bid_msgarch_var_025.csv")
write.csv(bid_msgarch_es_025,"bid_msgarch_es_025.csv")

msgarch_vol= read.csv("msgarch_vol.csv")
msgarch_vol= msgarch_vol[,-1]

bid_msgarch_var_01= read.csv("bid_msgarch_var_01.csv")
bid_msgarch_var_01= bid_msgarch_var_01[,-1]
bid_msgarch_es_01= read.csv("bid_msgarch_es_01.csv")
bid_msgarch_es_01= bid_msgarch_es_01[,-1]

bid_msgarch_var_05= read.csv("bid_msgarch_var_05.csv")
bid_msgarch_var_05= bid_msgarch_var_05[,-1]
bid_msgarch_es_05= read.csv("bid_msgarch_es_05.csv")
bid_msgarch_es_05= bid_msgarch_es_05[,-1]

bid_msgarch_var_025= read.csv("bid_msgarch_var_025.csv")
bid_msgarch_var_025= bid_msgarch_var_025[,-1]
bid_msgarch_es_025= read.csv("bid_msgarch_es_025.csv")
bid_msgarch_es_025= bid_msgarch_es_025[,-1]





# In sample estimate For C-RV with normal copula and the scaling factor estimation

m1= 1000
insam_est= c(rep(0,m1))

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
  insam_est[q]= mean(my_ecdf_inv_1)
}

insam_est= insam_est
write.csv(insam_est,"insam_est.csv")
insam_est= read.csv("insam_est.csv")
insam_est= insam_est[,-1]
scaling_factor= var(log_ret[1:1000])/mean(insam_est)



# Rolling window forecast for C-RV Model with normal copula

nor_mean_1= c(rep(0,m))

# benchmark model

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

bid_vol_nor= (nor_mean_1)
write.csv(bid_vol_nor,"bid_vol_nor.csv")
bid_vol_nor= read.csv("bid_vol_nor.csv")
bid_vol_nor= bid_vol_nor[,-1]
adj_bid_vol_nor= scaling_factor*(bid_vol_nor)
adj_bid_vol_nor= sqrt(adj_bid_vol_nor)




# In sample estimate For C-RV with t copula and the scaling factor estimation

m1= 1000
insam_est_t= c(rep(0,m1))

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
  insam_est_t[q]= mean(my_ecdf_inv_t)
}

insam_est_t= insam_est_t
write.csv(insam_est_t,"insam_est_t.csv")
insam_est_t= read.csv("insam_est_t.csv")
insam_est_t= insam_est_t[,-1]
scaling_factor_t= var(log_ret[1:1000])/mean(insam_est_t)


# Rolling window forecast for C-RV Model with t copula

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

bid_t_vol= (t_mean_1)
write.csv(bid_t_vol,"bid_t_vol.csv")
bid_t_vol= read.csv("bid_t_vol.csv")
bid_t_vol= bid_t_vol[,-1]
adj_bid_vol_t= (scaling_factor_t*bid_t_vol)
adj_bid_vol_t= sqrt(adj_bid_vol_t)


# in sample estimate CARR model and scaling factor

in_sam_ran= range[1:1000]
spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
e_garch_range <- ugarchfit(in_sam_ran, spec=spec.gjrGARCH_range)
est_range= e_garch_range@fit$var
scaling_carr= sd(log_ret[1:1000])/mean(est_range)


# Rolling window forecast by CARR model
carr_vol= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  train_fast1= range[start:end]
  spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
  e_garch_range <- ugarchfit(train_fast1, spec=spec.gjrGARCH_range)
  predi_garch_range= ugarchforecast(e_garch_range,n.ahead=1,data=train_fast1)
  fore= (predi_garch_range@forecast)
  carr_vol[q]= as.numeric(fore$sigma)
}

bid_carr_vol= (carr_vol)
write.csv(bid_carr_vol,"bid_carr_vol.csv")
bid_carr_vol= read.csv("bid_carr_vol.csv")
bid_carr_vol= bid_carr_vol[,-1]
bid_carr_vol= (bid_carr_vol)^2
adj_bid_carr_vol= scaling_carr*(bid_carr_vol)

#  Rolling window forecast for RobGARCHBoot Model

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

bid_rob_garch_var_01= (VaR_rob_01)
bid_rob_garch_var_05= (VaR_rob_05)
bid_rob_garch_var_025= (VaR_rob_025)

bid_rob_garch_ES_01= (ES_rob_01)
bid_rob_garch_ES_05= (ES_rob_05)
bid_rob_garch_ES_025= (ES_rob_025)

write.csv(bid_rob_garch_var_01,"bid_rob_garch_var_01.csv")
write.csv(bid_rob_garch_var_05,"bid_rob_garch_var_05.csv")
write.csv(bid_rob_garch_var_025,"bid_rob_garch_var_025.csv")

write.csv(bid_rob_garch_ES_01,"bid_rob_garch_ES_01.csv")
write.csv(bid_rob_garch_ES_05,"bid_rob_garch_ES_05.csv")
write.csv(bid_rob_garch_ES_025,"bid_rob_garch_ES_025.csv")



bid_rob_garch_var_01= read.csv("bid_rob_garch_var_01.csv")
bid_rob_garch_var_01= bid_rob_garch_var_01[,-1]
bid_rob_garch_var_025= read.csv("bid_rob_garch_var_025.csv")
bid_rob_garch_var_025= bid_rob_garch_var_025[,-1]
bid_rob_garch_var_05= read.csv("bid_rob_garch_var_05.csv")
bid_rob_garch_var_05= bid_rob_garch_var_05[,-1]

bid_rob_garch_es_01= read.csv("bid_rob_garch_es_01.csv")
bid_rob_garch_es_01= bid_rob_garch_es_01[,-1]
bid_rob_garch_es_025= read.csv("bid_rob_garch_es_025.csv")
bid_rob_garch_es_025= bid_rob_garch_es_025[,-1]
bid_rob_garch_es_05= read.csv("bid_rob_garch_es_05.csv")
bid_rob_garch_es_05= bid_rob_garch_es_05[,-1]


# GAS Model

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

bid_gas_var_01=  gas_var_01
bid_gas_es_01=  gas_es_01
bid_gas_var_05=  gas_var_05
bid_gas_es_05=  gas_es_05
bid_gas_var_025=  gas_var_025
bid_gas_es_025=  gas_es_025

write.csv(bid_gas_var_01,"bid_gas_var_01.csv")
write.csv(bid_gas_es_01,"bid_gas_es_01.csv")
write.csv(bid_gas_var_05,"bid_gas_var_05.csv")
write.csv(bid_gas_es_05,"bid_gas_es_05.csv")
write.csv(bid_gas_var_025,"bid_gas_var_025.csv")
write.csv(bid_gas_es_025,"bid_gas_es_025.csv")

bid_gas_var_01= read.csv("bid_gas_var_01.csv")
bid_gas_var_01= bid_gas_var_01[,-1]
bid_gas_es_01= read.csv("bid_gas_es_01.csv")
bid_gas_es_01= bid_gas_es_01[,-1]
bid_gas_var_05= read.csv("bid_gas_var_05.csv")
bid_gas_var_05= bid_gas_var_05[,-1]
bid_gas_es_05= read.csv("bid_gas_es_05.csv")
bid_gas_es_05= bid_gas_es_05[,-1]
bid_gas_var_025= read.csv("bid_gas_var_025.csv")
bid_gas_var_025= bid_gas_var_025[,-1]
bid_gas_es_025= read.csv("bid_gas_es_025.csv")
bid_gas_es_025= bid_gas_es_025[,-1]

# GARCH-ged Model

pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

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

garch_mean_ged= fit_mean_garch
bid_garch_ged= pre_garch_fast1
write.csv(bid_garch_ged,"bid_garch_ged.csv")
write.csv(garch_mean_ged,"garch_mean_ged.csv")
bid_garch_ged= read.csv("bid_garch_ged.csv")
bid_garch_ged= bid_garch_ged[,-1]
bid_mean_ged= read.csv("garch_mean_ged.csv")
bid_mean_ged= bid_mean_ged[,-1]

# Rolling window forecast for NAGARCH-ged Model


pre_nagarch_fast1= c(rep(0,m))

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
  pre_nagarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

bid_nagarch_ged= pre_nagarch_fast1
write.csv(bid_nagarch_ged,"bid_nagarch_ged.csv")
bid_nagarch_ged= read.csv("bid_nagarch_ged.csv")
bid_nagarch_ged= bid_nagarch_ged[,-1]

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

bid_gjrgarch_ged= pre_gjrgarch_fast1
write.csv(bid_gjrgarch_ged,"bid_gjrgarch_ged.csv")
bid_gjrgarch_ged= read.csv("bid_gjrgarch_ged.csv")
bid_gjrgarch_ged= bid_gjrgarch_ged[,-1]


# GARCH-sged Model

pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  fit_mean_garch[q]= as.numeric(pre_mean_garch$pred)
  spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sged")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH)
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_garch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

garch_mean_sged= fit_mean_garch
bid_garch_sged= pre_garch_fast1
write.csv(bid_garch_sged,"bid_garch_sged.csv")
write.csv(garch_mean_sged,"garch_mean_sged.csv")
bid_garch_sged= read.csv("bid_garch_sged.csv")
bid_garch_sged= bid_garch_sged[,-1]
bid_mean_sged= read.csv("garch_mean_sged.csv")
bid_mean_sged= bid_mean_sged[,-1]


# Rolling window forecast for NAGARCH-sged Model


pre_nagarch_fast1= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  
  spec.gjrGARCH = ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(1,1),submodel= "NAGARCH"), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sged")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH,solver = "hybrid")
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_nagarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

bid_nagarch_sged= pre_nagarch_fast1
write.csv(bid_nagarch_sged,"bid_nagarch_sged.csv")
bid_nagarch_sged= read.csv("bid_nagarch_sged.csv")
bid_nagarch_sged= bid_nagarch_sged[,-1]


#GJR-GARCH-sged Model

pre_gjrgarch_fast1=c(rep(0,m))
for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  res_gjrgarch= fit$residuals
  spec.garch = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sged")
  e_gjrgarch <- ugarchfit(res_gjrgarch, spec=spec.garch)
  e_forecast= ugarchforecast(e_gjrgarch,n.ahead=1,data=res_gjrgarch)
  fore= e_forecast@forecast
  pre_gjrgarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

bid_gjrgarch_sged= pre_gjrgarch_fast1
write.csv(bid_gjrgarch_sged,"bid_gjrgarch_sged.csv")
bid_gjrgarch_sged= read.csv("bid_gjrgarch_sged.csv")
bid_gjrgarch_sged= bid_gjrgarch_sged[,-1]


# GARCH-jsu Model

pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  fit_mean_garch[q]= as.numeric(pre_mean_garch$pred)
  spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="jsu")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH)
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_garch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

garch_mean_jsu= fit_mean_garch
bid_garch_jsu= pre_garch_fast1
write.csv(bid_garch_jsu,"bid_garch_jsu.csv")
write.csv(garch_mean_jsu,"garch_mean_jsu.csv")
bid_garch_jsu= read.csv("bid_garch_jsu.csv")
bid_garch_jsu= bid_garch_jsu[,-1]
bid_mean_jsu= read.csv("garch_mean_jsu.csv")
bid_mean_jsu= bid_mean_jsu[,-1]


# Rolling window forecast for NAGARCH-jsu Model


pre_nagarch_fast1= c(rep(0,m))

for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  
  spec.gjrGARCH = ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(1,1),submodel= "NAGARCH"), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="jsu")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH,solver = "hybrid")
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_nagarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

bid_nagarch_jsu= pre_nagarch_fast1
write.csv(bid_nagarch_jsu,"bid_nagarch_jsu.csv")
bid_nagarch_jsu= read.csv("bid_nagarch_jsu.csv")
bid_nagarch_jsu= bid_nagarch_jsu[,-1]


#GJR-GARCH-jsu Model

pre_gjrgarch_fast1=c(rep(0,m))
for(q in 1:m){
  
  start= q
  end= 999+q
  
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 3,max.q = 0)
  res_gjrgarch= fit$residuals
  spec.garch = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="jsu")
  e_gjrgarch <- ugarchfit(res_gjrgarch, spec=spec.garch)
  e_forecast= ugarchforecast(e_gjrgarch,n.ahead=1,data=res_gjrgarch)
  fore= e_forecast@forecast
  pre_gjrgarch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

bid_gjrgarch_jsu= pre_gjrgarch_fast1
write.csv(bid_gjrgarch_jsu,"bid_gjrgarch_jsu.csv")
bid_gjrgarch_jsu= read.csv("bid_gjrgarch_jsu.csv")
bid_gjrgarch_jsu= bid_gjrgarch_jsu[,-1]


# VaR assuming sstd distribution on the innovations

fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 3,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="sstd")
e_garch_sstd <- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)


# VaR_sstd for GARCH model
VaR_sstd_garch =  (bid_garch)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                    shape=coef(e_garch_sstd)["shape"])

# VaR_sstd for NAGARCH model

VaR_sstd_nagarch =  (bid_nagarch)*qdist("sstd", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                                        shape=coef(e_garch_sstd)["shape"])

# VaR_sstd for C-RV with normal copula model


VaR_sstd_nor_cop =  (adj_bid_vol_nor)*qdist("sstd", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                                            shape=coef(e_garch_sstd)["shape"])

# VaR_sstd for C-RV with t copula model

VaR_sstd_t_cop =  (adj_bid_vol_t)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                        shape=coef(e_garch_sstd)["shape"])


# VaR_sstd for CARR model

VaR_sstd_carr =  (adj_bid_carr_vol)*qdist("sstd", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                                          shape=coef(e_garch_sstd)["shape"])
# VaR_sstd GJR-GARCH model
VaR_sstd_gjrgarch =  (bid_gjrgarch)*qdist("sstd", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sstd)["skew"],
                                          shape=coef(e_garch_sstd)["shape"])

# VaR with ged on the innovations
fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 3,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="ged")
e_garch <- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)


# VaR_ged for GARCH model
VaR_ged_garch =  (bid_garch_ged)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                       shape=coef(e_garch)["shape"])

# VaR_ged for NAGARCH model

VaR_ged_nagarch =  (bid_nagarch_ged)*qdist("ged", p=0.05, mu = 0, sigma = 1,
                                           shape=coef(e_garch)["shape"])

# VaR_ged for C-RV with normal copula model


VaR_ged_nor_cop =  (adj_bid_vol_nor)*qdist("ged", p=0.05, mu = 0, sigma = 1,
                                           shape=coef(e_garch)["shape"])

# VaR_ged for C-RV with t copula model

VaR_ged_t_cop =  (adj_bid_vol_t)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                       shape=coef(e_garch)["shape"])


# VaR_ged for CARR model

VaR_ged_carr =  (adj_bid_carr_vol)*qdist("ged", p=0.05, mu = 0, sigma = 1,
                                         shape=coef(e_garch)["shape"])
# VaR_ged GJR-GARCH model
VaR_ged_gjrgarch =  (bid_gjrgarch_ged)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                             shape=coef(e_garch)["shape"])


# VaR assuming sged distribution on the innovations

fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 3,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="sged")
e_garch_sged <- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)

# VaR_sged for GARCH model
VaR_sged_garch =  (bid_garch_sged)*qdist("sged", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sged)["skew"],
                                         shape=coef(e_garch_sged)["shape"])

# VaR_sged for NAGARCH model

VaR_sged_nagarch =  (bid_nagarch_sged)*qdist("sged", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                                             shape=coef(e_garch_sged)["shape"])

# VaR_sged for C-RV with normal copula model


VaR_sged_nor_cop =  (adj_bid_vol_nor)*qdist("sged", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                                            shape=coef(e_garch_sged)["shape"])

# VaR_sged for C-RV with t copula model

VaR_sged_t_cop =  (adj_bid_vol_t)*qdist("sged", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sged)["skew"],
                                        shape=coef(e_garch_sged)["shape"])


# VaR_sged for CARR model

VaR_sged_carr =  (adj_bid_carr_vol)*qdist("sged", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                                          shape=coef(e_garch_sged)["shape"])
# VaR_sged GJR-GARCH model
VaR_sged_gjrgarch =  (bid_gjrgarch_sged)*qdist("sged", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_sged)["skew"],
                                               shape=coef(e_garch_sged)["shape"])


# VaR assuming jsu distribution on the innovations

fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 3,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="jsu")
e_garch_jsu <- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)

# VaR_jsu for GARCH model
VaR_jsu_garch =  (bid_garch_jsu)*qdist("jsu", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_jsu)["skew"],
                                       shape=coef(e_garch_jsu)["shape"])

# VaR_jsu for NAGARCH model

VaR_jsu_nagarch =  (bid_nagarch_jsu)*qdist("jsu", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                                           shape=coef(e_garch_jsu)["shape"])

# VaR_jsu for C-RV with normal copula model


VaR_jsu_nor_cop =  (adj_bid_vol_nor)*qdist("jsu", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                                           shape=coef(e_garch_jsu)["shape"])

# VaR_jsu for C-RV with t copula model

VaR_jsu_t_cop =  (adj_bid_vol_t)*qdist("jsu", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_jsu)["skew"],
                                       shape=coef(e_garch_jsu)["shape"])


# VaR_jsu for CARR model

VaR_jsu_carr =  (adj_bid_carr_vol)*qdist("jsu", p=0.05, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                                         shape=coef(e_garch_jsu)["shape"])
# VaR_jsu GJR-GARCH model
VaR_jsu_gjrgarch =  (bid_gjrgarch_jsu)*qdist("jsu", p=0.05, mu = 0, sigma = 1, skew=coef(e_garch_jsu)["skew"],
                                             shape=coef(e_garch_jsu)["shape"])



test= log_ret[1001:2343]
actual= test-bid_mean

# VaR Backtesting

backtest_sstd_nor_cop= BacktestVaR(data = actual, VaR = VaR_sstd_nor_cop, alpha = 0.05)

backtest_sstd_garch= BacktestVaR(data = actual, VaR = VaR_sstd_garch, alpha = 0.05)

backtest_sstd_carr= BacktestVaR(data = actual, VaR = VaR_sstd_carr, alpha = 0.05)

backtest_sstd_t_cop= BacktestVaR(data = actual, VaR = VaR_sstd_t_cop, alpha = 0.05)

backtest_robgarch= BacktestVaR(data = actual, VaR = bid_rob_garch_var_05, alpha = 0.05)

backtest_sstd_nagarch= BacktestVaR(data = actual, VaR = VaR_sstd_nagarch, alpha = 0.05)

backtest_ged_nor_cop= BacktestVaR(data = actual, VaR = VaR_ged_nor_cop, alpha = 0.05)

backtest_ged_garch= BacktestVaR(data = actual, VaR = VaR_ged_garch, alpha = 0.05)

backtest_ged_carr= BacktestVaR(data = actual, VaR = VaR_ged_carr, alpha = 0.05)

backtest_ged_t_cop= BacktestVaR(data = actual, VaR = VaR_ged_t_cop, alpha = 0.05)
backtest_ged_nagarch= BacktestVaR(data = actual, VaR = VaR_ged_nagarch, alpha = 0.05)
backtest_ged_gjrgarch= BacktestVaR(data = actual, VaR = VaR_ged_gjrgarch, alpha = 0.05)

backtest_msgarch= BacktestVaR(data = actual, VaR = bid_msgarch_var_05, alpha = 0.05)

backtest_gas= BacktestVaR(data = actual, VaR = bid_gas_var_05, alpha = 0.05)

backtest_sstd_gjrgarch= BacktestVaR(data = actual, VaR = VaR_sstd_gjrgarch, alpha = 0.05)

backtest_sged_nor_cop= BacktestVaR(data = actual, VaR = VaR_sged_nor_cop, alpha = 0.05)

backtest_sged_garch= BacktestVaR(data = actual, VaR = VaR_sged_garch, alpha = 0.05)

backtest_sged_carr= BacktestVaR(data = actual, VaR = VaR_sged_carr, alpha = 0.05)

backtest_sged_t_cop= BacktestVaR(data = actual, VaR = VaR_sged_t_cop, alpha = 0.05)
backtest_sged_nagarch= BacktestVaR(data = actual, VaR = VaR_sged_nagarch, alpha = 0.05)
backtest_sged_gjrgarch= BacktestVaR(data = actual, VaR = VaR_sged_gjrgarch, alpha = 0.05)

backtest_jsu_nor_cop= BacktestVaR(data = actual, VaR = VaR_jsu_nor_cop, alpha = 0.05)

backtest_jsu_garch= BacktestVaR(data = actual, VaR = VaR_jsu_garch, alpha = 0.05)

backtest_jsu_carr= BacktestVaR(data = actual, VaR = VaR_jsu_carr, alpha = 0.05)

backtest_jsu_t_cop= BacktestVaR(data = actual, VaR = VaR_jsu_t_cop, alpha = 0.05)
backtest_jsu_nagarch= BacktestVaR(data = actual, VaR = VaR_jsu_nagarch, alpha = 0.05)
backtest_jsu_gjrgarch= BacktestVaR(data = actual, VaR = VaR_jsu_gjrgarch, alpha = 0.05)


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

uc_sged_nor_cop= (backtest_sged_nor_cop$LRuc["Pvalue" ])
uc_sged_t_cop= (backtest_sged_t_cop$LRuc["Pvalue" ])
uc_sged_garch= backtest_sged_garch$LRuc["Pvalue" ]
uc_sged_nagarch= backtest_sged_nagarch$LRuc["Pvalue" ]
uc_sged_gjrgarch= backtest_sged_gjrgarch$LRuc["Pvalue" ]
uc_sged_carr= backtest_sged_carr$LRuc["Pvalue" ]
uc_jsu_nor_cop= (backtest_jsu_nor_cop$LRuc["Pvalue" ])
uc_jsu_t_cop= (backtest_jsu_t_cop$LRuc["Pvalue" ])
uc_jsu_garch= backtest_jsu_garch$LRuc["Pvalue" ]
uc_jsu_nagarch= backtest_jsu_nagarch$LRuc["Pvalue" ]
uc_jsu_gjrgarch= backtest_jsu_gjrgarch$LRuc["Pvalue" ]
uc_jsu_carr= backtest_jsu_carr$LRuc["Pvalue" ]


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

cc_sged_nor_cop= (backtest_sged_nor_cop$LRcc["Pvalue" ])
cc_sged_t_cop= (backtest_sged_t_cop$LRcc["Pvalue" ])
cc_sged_garch= backtest_sged_garch$LRcc["Pvalue" ]
cc_sged_nagarch= backtest_sged_nagarch$LRcc["Pvalue" ]
cc_sged_gjrgarch= backtest_sged_gjrgarch$LRcc["Pvalue" ]
cc_sged_carr= backtest_sged_carr$LRcc["Pvalue" ]

cc_jsu_nor_cop= (backtest_jsu_nor_cop$LRcc["Pvalue" ])
cc_jsu_t_cop= (backtest_jsu_t_cop$LRcc["Pvalue" ])
cc_jsu_garch= backtest_jsu_garch$LRcc["Pvalue" ]
cc_jsu_nagarch= backtest_jsu_nagarch$LRcc["Pvalue" ]
cc_jsu_gjrgarch= backtest_jsu_gjrgarch$LRcc["Pvalue" ]
cc_jsu_carr= backtest_jsu_carr$LRcc["Pvalue" ]



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

DQ_sged_nor_cop= as.numeric(backtest_sged_nor_cop$DQ$pvalue)
DQ_sged_t_cop= as.numeric(backtest_sged_t_cop$DQ$pvalue)
DQ_sged_garch=as.numeric( backtest_sged_garch$DQ$pvalue)
DQ_sged_nagarch=as.numeric( backtest_sged_nagarch$DQ$pvalue)
DQ_sged_gjrgarch=as.numeric( backtest_sged_gjrgarch$DQ$pvalue)
DQ_sged_carr=as.numeric( backtest_sged_carr$DQ$pvalue)

DQ_jsu_nor_cop= as.numeric(backtest_jsu_nor_cop$DQ$pvalue)
DQ_jsu_t_cop= as.numeric(backtest_jsu_t_cop$DQ$pvalue)
DQ_jsu_garch=as.numeric( backtest_jsu_garch$DQ$pvalue)
DQ_jsu_nagarch=as.numeric( backtest_jsu_nagarch$DQ$pvalue)
DQ_jsu_gjrgarch=as.numeric( backtest_jsu_gjrgarch$DQ$pvalue)
DQ_jsu_carr=as.numeric( backtest_jsu_carr$DQ$pvalue)

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

QL_sged_nor_cop= as.numeric(backtest_sged_nor_cop$Loss$Loss)
QL_sged_t_cop= as.numeric(backtest_sged_t_cop$Loss$Loss)
QL_sged_garch=as.numeric( backtest_sged_garch$Loss$Loss)
QL_sged_nagarch=as.numeric( backtest_sged_nagarch$Loss$Loss)
QL_sged_gjrgarch=as.numeric( backtest_sged_gjrgarch$Loss$Loss)
QL_sged_carr=as.numeric( backtest_sged_carr$Loss$Loss)

QL_jsu_nor_cop= as.numeric(backtest_jsu_nor_cop$Loss$Loss)
QL_jsu_t_cop= as.numeric(backtest_jsu_t_cop$Loss$Loss)
QL_jsu_garch=as.numeric( backtest_jsu_garch$Loss$Loss)
QL_jsu_nagarch=as.numeric( backtest_jsu_nagarch$Loss$Loss)
QL_jsu_gjrgarch=as.numeric( backtest_jsu_gjrgarch$Loss$Loss)
QL_jsu_carr=as.numeric( backtest_jsu_carr$Loss$Loss)

# AE Test
AE_ged_nor_cop= as.numeric(backtest_ged_nor_cop$AE)
AE_ged_t_cop= as.numeric(backtest_ged_t_cop$AE)
AE_ged_garch=as.numeric( backtest_ged_garch$AE)
AE_ged_nagarch=as.numeric( backtest_ged_nagarch$AE)
AE_ged_gjrgarch=as.numeric( backtest_ged_gjrgarch$AE)
AE_ged_carr=as.numeric( backtest_ged_carr$AE)
AE_msgarch=as.numeric( backtest_msgarch$AE)
AE_robgarch=as.numeric( backtest_robgarch$AE)
AE_gas=as.numeric( backtest_gas$AE)
AE_sstd_nor_cop= as.numeric(backtest_sstd_nor_cop$AE)
AE_sstd_t_cop= as.numeric(backtest_sstd_t_cop$AE)
AE_sstd_garch=as.numeric( backtest_sstd_garch$AE)
AE_sstd_nagarch=as.numeric( backtest_sstd_nagarch$AE)
AE_sstd_gjrgarch=as.numeric( backtest_sstd_gjrgarch$AE)
AE_sstd_carr=as.numeric( backtest_sstd_carr$AE)

AE_sged_nor_cop= as.numeric(backtest_sged_nor_cop$AE)
AE_sged_t_cop= as.numeric(backtest_sged_t_cop$AE)
AE_sged_garch=as.numeric( backtest_sged_garch$AE)
AE_sged_nagarch=as.numeric( backtest_sged_nagarch$AE)
AE_sged_gjrgarch=as.numeric( backtest_sged_gjrgarch$AE)
AE_sged_carr=as.numeric( backtest_sged_carr$AE)

AE_jsu_nor_cop= as.numeric(backtest_jsu_nor_cop$AE)
AE_jsu_t_cop= as.numeric(backtest_jsu_t_cop$AE)
AE_jsu_garch=as.numeric( backtest_jsu_garch$AE)
AE_jsu_nagarch=as.numeric( backtest_jsu_nagarch$AE)
AE_jsu_gjrgarch=as.numeric( backtest_jsu_gjrgarch$AE)
AE_jsu_carr=as.numeric( backtest_jsu_carr$AE)

# Expected Shortfall and McF Test

f1 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])

ES_sstd_garch = (bid_garch)*integrate(f1, 0, 0.05)$value/0.05

# McF test
mcf_sstd_garch= (ESTest(0.05, actual, ES_sstd_garch, VaR_sstd_garch, boot = TRUE))



f2 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_nor_vol = (adj_bid_vol_nor)*integrate(f2, 0, 0.05)$value/0.05
# McF test
mcf_sstd_nor_cop= (ESTest(0.05, actual, ES_sstd_nor_vol, VaR_sstd_nor_cop, boot = TRUE))



f4 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_t_vol = (adj_bid_vol_t)*integrate(f4, 0, 0.05)$value/0.05
# McF test

mcf_sstd_t_cop=(ESTest(0.05, actual, ES_sstd_t_vol,VaR_sstd_t_cop, boot = TRUE))

# McF test

mcf_robgarch =(ESTest(0.05, actual, bid_rob_garch_es_05, bid_rob_garch_var_05, boot = TRUE))


f3 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_carr = (adj_bid_carr_vol)*integrate(f3, 0, 0.05)$value/0.05
# McF test
mcf_sstd_carr= (ESTest(0.05, actual, ES_sstd_carr, VaR_sstd_carr, boot = TRUE))

f5 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])
ES_sstd_nagarch = (bid_nagarch)*integrate(f5, 0, 0.05)$value/0.05
# McF test
mcf_sstd_nagarch= (ESTest(0.05, actual, ES_sstd_nagarch, VaR_sstd_nagarch, boot = TRUE))

# McF test
mcf_msgarch=(ESTest(0.05, actual, bid_msgarch_es_05, bid_msgarch_var_05, boot = TRUE))

mcf_gas =(ESTest(0.05, test, bid_gas_es_05, bid_gas_var_05, boot = TRUE))

f6 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sstd)["skew"],
                       shape=coef(e_garch_sstd)["shape"])

ES_sstd_gjrgarch = (bid_gjrgarch)*integrate(f6, 0, 0.05)$value/0.05

# McF test
mcf_sstd_gjrgarch=(ESTest(0.05, actual, ES_sstd_gjrgarch, VaR_sstd_gjrgarch, boot = TRUE))

f1_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])

ES_ged_garch = (bid_garch_ged)*integrate(f1_ged, 0, 0.05)$value/0.05

# McF test
mcf_ged_garch= (ESTest(0.05, actual, ES_ged_garch, VaR_ged_garch, boot = TRUE))



f2_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_nor_vol = (adj_bid_vol_nor)*integrate(f2_ged, 0, 0.05)$value/0.05
# McF test
mcf_ged_nor_cop= (ESTest(0.05, actual, ES_ged_nor_vol, VaR_ged_nor_cop, boot = TRUE))



f4_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_t_vol = (adj_bid_vol_t)*integrate(f4_ged, 0, 0.05)$value/0.05
# McF test

mcf_ged_t_cop=(ESTest(0.05, actual, ES_ged_t_vol,VaR_ged_t_cop, boot = TRUE))

f3_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_carr = (adj_bid_carr_vol)*integrate(f3_ged, 0, 0.05)$value/0.05
# McF test
mcf_ged_carr= (ESTest(0.05, actual, ES_ged_carr, VaR_ged_carr, boot = TRUE))

f5_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])
ES_ged_nagarch = (bid_nagarch_ged)*integrate(f5_ged, 0, 0.05)$value/0.05
# McF test
mcf_ged_nagarch= (ESTest(0.05, actual, ES_ged_nagarch, VaR_ged_nagarch, boot = TRUE))

f6_ged = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                           shape=coef(e_garch)["shape"])

ES_ged_gjrgarch = (bid_gjrgarch_ged)*integrate(f6_ged, 0, 0.05)$value/0.05

# McF test
mcf_ged_gjrgarch=(ESTest(0.05, actual, ES_ged_gjrgarch, VaR_ged_gjrgarch, boot = TRUE))

# Expected Shortfall and McF Test

f1 = function(x) qdist("sged", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                       shape=coef(e_garch_sged)["shape"])

ES_sged_garch = (bid_garch_sged)*integrate(f1, 0, 0.05)$value/0.05

# McF test
mcf_sged_garch= (ESTest(0.05, actual, ES_sged_garch, VaR_sged_garch, boot = TRUE))



f2 = function(x) qdist("sged", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                       shape=coef(e_garch_sged)["shape"])
ES_sged_nor_vol = (adj_bid_vol_nor)*integrate(f2, 0, 0.05)$value/0.05
# McF test
mcf_sged_nor_cop= (ESTest(0.05, actual, ES_sged_nor_vol, VaR_sged_nor_cop, boot = TRUE))



f4 = function(x) qdist("sged", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                       shape=coef(e_garch_sged)["shape"])
ES_sged_t_vol = (adj_bid_vol_t)*integrate(f4, 0, 0.05)$value/0.05
# McF test

mcf_sged_t_cop=(ESTest(0.05, actual, ES_sged_t_vol,VaR_sged_t_cop, boot = TRUE))

f3 = function(x) qdist("sged", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                       shape=coef(e_garch_sged)["shape"])
ES_sged_carr = (adj_bid_carr_vol)*integrate(f3, 0, 0.05)$value/0.05
# McF test
mcf_sged_carr= (ESTest(0.05, actual, ES_sged_carr, VaR_sged_carr, boot = TRUE))

f5 = function(x) qdist("sged", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                       shape=coef(e_garch_sged)["shape"])
ES_sged_nagarch = (bid_nagarch_sged)*integrate(f5, 0, 0.05)$value/0.05
# McF test
mcf_sged_nagarch= (ESTest(0.05, actual, ES_sged_nagarch, VaR_sged_nagarch, boot = TRUE))

f6 = function(x) qdist("sged", p=x, mu = 0, sigma = 1,skew=coef(e_garch_sged)["skew"],
                       shape=coef(e_garch_sged)["shape"])

ES_sged_gjrgarch = (bid_gjrgarch_sged)*integrate(f6, 0, 0.05)$value/0.05

# McF test
mcf_sged_gjrgarch=(ESTest(0.05, actual, ES_sged_gjrgarch, VaR_sged_gjrgarch, boot = TRUE))
# Expected Shortfall and McF Test

f1 = function(x) qdist("jsu", p=x, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                       shape=coef(e_garch_jsu)["shape"])

ES_jsu_garch = (bid_garch_jsu)*integrate(f1, 0, 0.05)$value/0.05

# McF test
mcf_jsu_garch= (ESTest(0.05, actual, ES_jsu_garch, VaR_jsu_garch, boot = TRUE))



f2 = function(x) qdist("jsu", p=x, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                       shape=coef(e_garch_jsu)["shape"])
ES_jsu_nor_vol = (adj_bid_vol_nor)*integrate(f2, 0, 0.05)$value/0.05
# McF test
mcf_jsu_nor_cop= (ESTest(0.05, actual, ES_jsu_nor_vol, VaR_jsu_nor_cop, boot = TRUE))



f4 = function(x) qdist("jsu", p=x, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                       shape=coef(e_garch_jsu)["shape"])
ES_jsu_t_vol = (adj_bid_vol_t)*integrate(f4, 0, 0.05)$value/0.05
# McF test

mcf_jsu_t_cop=(ESTest(0.05, actual, ES_jsu_t_vol,VaR_jsu_t_cop, boot = TRUE))

f3 = function(x) qdist("jsu", p=x, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                       shape=coef(e_garch_jsu)["shape"])
ES_jsu_carr = (adj_bid_carr_vol)*integrate(f3, 0, 0.05)$value/0.05
# McF test
mcf_jsu_carr= (ESTest(0.05, actual, ES_jsu_carr, VaR_jsu_carr, boot = TRUE))

f5 = function(x) qdist("jsu", p=x, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                       shape=coef(e_garch_jsu)["shape"])
ES_jsu_nagarch = (bid_nagarch_jsu)*integrate(f5, 0, 0.05)$value/0.05
# McF test
mcf_jsu_nagarch= (ESTest(0.05, actual, ES_jsu_nagarch, VaR_jsu_nagarch, boot = TRUE))

f6 = function(x) qdist("jsu", p=x, mu = 0, sigma = 1,skew=coef(e_garch_jsu)["skew"],
                       shape=coef(e_garch_jsu)["shape"])

ES_jsu_gjrgarch = (bid_gjrgarch_jsu)*integrate(f6, 0, 0.05)$value/0.05

# McF test
mcf_jsu_gjrgarch=(ESTest(0.05, actual, ES_jsu_gjrgarch, VaR_jsu_gjrgarch, boot = TRUE))

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

hits_nor_sged_cop= (mcf_sged_nor_cop$actual.exceed*100)/1343
hits_t_sged_cop= (mcf_sged_t_cop$actual.exceed*100)/1343
hits_sged_garch= (mcf_sged_garch$actual.exceed*100)/1343
hits_sged_nagarch= (mcf_sged_nagarch$actual.exceed*100)/1343
hits_sged_gjrgarch= (mcf_sged_gjrgarch$actual.exceed*100)/1343
hits_sged_carr= (mcf_sged_carr$actual.exceed*100)/1343
hits_nor_jsu_cop= (mcf_jsu_nor_cop$actual.exceed*100)/1343
hits_t_jsu_cop= (mcf_jsu_t_cop$actual.exceed*100)/1343
hits_jsu_garch= (mcf_jsu_garch$actual.exceed*100)/1343
hits_jsu_nagarch= (mcf_jsu_nagarch$actual.exceed*100)/1343
hits_jsu_gjrgarch= (mcf_jsu_gjrgarch$actual.exceed*100)/1343
hits_jsu_carr= (mcf_jsu_carr$actual.exceed*100)/1343



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

McF_p_sged_nor_cop= mcf_sged_nor_cop$p.value
McF_p_sged_t_cop= mcf_sged_t_cop$p.value
McF_p_sged_garch= mcf_sged_garch$p.value
McF_p_sged_nagarch= mcf_sged_nagarch$p.value
McF_p_sged_gjrgarch= mcf_sged_gjrgarch$p.value
McF_p_sged_carr= mcf_sged_carr$p.value

McF_p_jsu_nor_cop= mcf_jsu_nor_cop$p.value
McF_p_jsu_t_cop= mcf_jsu_t_cop$p.value
McF_p_jsu_garch= mcf_jsu_garch$p.value
McF_p_jsu_nagarch= mcf_jsu_nagarch$p.value
McF_p_jsu_gjrgarch= mcf_jsu_gjrgarch$p.value
McF_p_jsu_carr= mcf_jsu_carr$p.value



# NF Test
NF_sstd_nor_cop= cc_backtest(actual, VaR_sstd_nor_cop, ES_sstd_nor_vol, alpha=0.05)
NF_sstd_t_cop= cc_backtest(actual, VaR_sstd_t_cop, ES_sstd_t_vol, alpha=0.05)
NF_sstd_garch= cc_backtest(actual, VaR_sstd_garch, ES_sstd_garch, alpha = 0.05)
NF_sstd_gjrgarch= cc_backtest(actual, VaR_sstd_gjrgarch, ES_sstd_gjrgarch,alpha = 0.05)
NF_sstd_nagarch= cc_backtest(actual, VaR_sstd_nagarch, ES_sstd_nagarch,alpha = 0.05)
NF_sstd_carr= cc_backtest(actual, VaR_sstd_carr, ES_sstd_carr, alpha=0.05)

NF_msgarch= cc_backtest(actual, bid_msgarch_var_05, bid_msgarch_es_05, alpha=0.05)
NF_gas= cc_backtest(actual, bid_gas_var_05, bid_gas_es_05, alpha=0.05)
NF_robgarch= cc_backtest(actual, bid_rob_garch_var_05, bid_rob_garch_es_05, alpha=0.05)
NF_ged_nor_cop= cc_backtest(actual, VaR_ged_nor_cop, ES_ged_nor_vol, alpha=0.05)
NF_ged_t_cop= cc_backtest(actual, VaR_ged_t_cop, ES_ged_t_vol, alpha=0.05)
NF_ged_garch= cc_backtest(actual, VaR_ged_garch, ES_ged_garch, alpha = 0.05)
NF_ged_gjrgarch= cc_backtest(actual, VaR_ged_gjrgarch, ES_ged_gjrgarch,alpha = 0.05)
NF_ged_nagarch= cc_backtest(actual, VaR_ged_nagarch, ES_ged_nagarch,alpha = 0.05)
NF_ged_carr= cc_backtest(actual, VaR_ged_carr, ES_ged_carr, alpha=0.05)

NF_sged_nor_cop= cc_backtest(actual, VaR_sged_nor_cop, ES_sged_nor_vol, alpha=0.05)
NF_sged_t_cop= cc_backtest(actual, VaR_sged_t_cop, ES_sged_t_vol, alpha=0.05)
NF_sged_garch= cc_backtest(actual, VaR_sged_garch, ES_sged_garch, alpha = 0.05)
NF_sged_gjrgarch= cc_backtest(actual, VaR_sged_gjrgarch, ES_sged_gjrgarch,alpha = 0.05)
NF_sged_nagarch= cc_backtest(actual, VaR_sged_nagarch, ES_sged_nagarch,alpha = 0.05)
NF_sged_carr= cc_backtest(actual, VaR_sged_carr, ES_sged_carr, alpha=0.05)
NF_jsu_nor_cop= cc_backtest(actual, VaR_jsu_nor_cop, ES_jsu_nor_vol, alpha=0.05)
NF_jsu_t_cop= cc_backtest(actual, VaR_jsu_t_cop, ES_jsu_t_vol, alpha=0.05)
NF_jsu_garch= cc_backtest(actual, VaR_jsu_garch, ES_jsu_garch, alpha = 0.05)
NF_jsu_gjrgarch= cc_backtest(actual, VaR_jsu_gjrgarch, ES_jsu_gjrgarch,alpha = 0.05)
NF_jsu_nagarch= cc_backtest(actual, VaR_jsu_nagarch, ES_jsu_nagarch,alpha = 0.05)
NF_jsu_carr= cc_backtest(actual, VaR_jsu_carr, ES_jsu_carr, alpha=0.05)

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

NF_p_sged_nor_cop= NF_sged_nor_cop$pvalue_twosided_simple
NF_p_sged_t_cop= NF_sged_t_cop$pvalue_twosided_simple
NF_p_sged_garch= NF_sged_garch$pvalue_twosided_simple
NF_p_sged_nagarch= NF_sged_nagarch$pvalue_twosided_simple
NF_p_sged_gjrgarch= NF_sged_gjrgarch$pvalue_twosided_simple
NF_p_sged_carr= NF_sged_carr$pvalue_twosided_simple

NF_p_jsu_nor_cop= NF_jsu_nor_cop$pvalue_twosided_simple
NF_p_jsu_t_cop= NF_jsu_t_cop$pvalue_twosided_simple
NF_p_jsu_garch= NF_jsu_garch$pvalue_twosided_simple
NF_p_jsu_nagarch= NF_jsu_nagarch$pvalue_twosided_simple
NF_p_jsu_gjrgarch= NF_jsu_gjrgarch$pvalue_twosided_simple
NF_p_jsu_carr= NF_jsu_carr$pvalue_twosided_simple

# BD Test

BD_sstd_nor_cop= esr_backtest(actual, VaR_sstd_nor_cop, ES_sstd_nor_vol, alpha=0.05,version = 1)
BD_sstd_t_cop= esr_backtest(actual, VaR_sstd_t_cop, ES_sstd_t_vol, alpha=0.05,version = 1)
BD_sstd_garch= esr_backtest(actual, VaR_sstd_garch, ES_sstd_garch, alpha = 0.05,version = 1)
BD_sstd_gjrgarch= esr_backtest(actual, VaR_sstd_gjrgarch, ES_sstd_gjrgarch,alpha = 0.05,version = 1)
BD_sstd_nagarch= esr_backtest(actual, VaR_sstd_nagarch, ES_sstd_nagarch,alpha = 0.05,version = 1)
BD_sstd_carr= esr_backtest(actual, VaR_sstd_carr, ES_sstd_carr, alpha=0.05,version = 1)
BD_msgarch= esr_backtest(actual, bid_msgarch_var_05, bid_msgarch_es_05, alpha=0.05,version = 1)
BD_gas= esr_backtest(actual, bid_gas_var_05, bid_gas_es_05, alpha=0.05,version = 1)
BD_robgarch= esr_backtest(actual, bid_rob_garch_var_05, bid_rob_garch_es_05, alpha=0.05,version = 1)
BD_ged_nor_cop= esr_backtest(actual, VaR_ged_nor_cop, ES_ged_nor_vol, alpha=0.05,version = 1)
BD_ged_t_cop= esr_backtest(actual, VaR_ged_t_cop, ES_ged_t_vol, alpha=0.05,version = 1)
BD_ged_garch= esr_backtest(actual, VaR_ged_garch, ES_ged_garch, alpha = 0.05,version = 1)
BD_ged_gjrgarch= esr_backtest(actual, VaR_ged_gjrgarch, ES_ged_gjrgarch,alpha = 0.05,version = 1)
BD_ged_nagarch= esr_backtest(actual, VaR_ged_nagarch, ES_ged_nagarch,alpha = 0.05,version = 1)
BD_ged_carr= esr_backtest(actual, VaR_ged_carr, ES_ged_carr, alpha=0.05,version = 1)

BD_sged_nor_cop= esr_backtest(actual, VaR_sged_nor_cop, ES_sged_nor_vol, alpha=0.05,version = 1)
BD_sged_t_cop= esr_backtest(actual, VaR_sged_t_cop, ES_sged_t_vol, alpha=0.05,version = 1)
BD_sged_garch= esr_backtest(actual, VaR_sged_garch, ES_sged_garch, alpha = 0.05,version = 1)
BD_sged_gjrgarch= esr_backtest(actual, VaR_sged_gjrgarch, ES_sged_gjrgarch,alpha = 0.05,version = 1)
BD_sged_nagarch= esr_backtest(actual, VaR_sged_nagarch, ES_sged_nagarch,alpha = 0.05,version = 1)
BD_sged_carr= esr_backtest(actual, VaR_sged_carr, ES_sged_carr, alpha=0.05,version = 1)

BD_jsu_nor_cop= esr_backtest(actual, VaR_jsu_nor_cop, ES_jsu_nor_vol, alpha=0.05,version = 1)
BD_jsu_t_cop= esr_backtest(actual, VaR_jsu_t_cop, ES_jsu_t_vol, alpha=0.05,version = 1)
BD_jsu_garch= esr_backtest(actual, VaR_jsu_garch, ES_jsu_garch, alpha = 0.05,version = 1)
BD_jsu_gjrgarch= esr_backtest(actual, VaR_jsu_gjrgarch, ES_jsu_gjrgarch,alpha = 0.05,version = 1)
BD_jsu_nagarch= esr_backtest(actual, VaR_jsu_nagarch, ES_jsu_nagarch,alpha = 0.05,version = 1)
BD_jsu_carr= esr_backtest(actual, VaR_jsu_carr, ES_jsu_carr, alpha=0.05,version = 1)

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

BD_p_sged_nor_cop= BD_sged_nor_cop$pvalue_twosided_asymptotic
BD_p_sged_t_cop= BD_sged_t_cop$pvalue_twosided_asymptotic
BD_p_sged_garch= BD_sged_garch$pvalue_twosided_asymptotic
BD_p_sged_nagarch= BD_sged_nagarch$pvalue_twosided_asymptotic
BD_p_sged_gjrgarch= BD_sged_gjrgarch$pvalue_twosided_asymptotic
BD_p_sged_carr= BD_sged_carr$pvalue_twosided_asymptotic

BD_p_jsu_nor_cop= BD_jsu_nor_cop$pvalue_twosided_asymptotic
BD_p_jsu_t_cop= BD_jsu_t_cop$pvalue_twosided_asymptotic
BD_p_jsu_garch= BD_jsu_garch$pvalue_twosided_asymptotic
BD_p_jsu_nagarch= BD_jsu_nagarch$pvalue_twosided_asymptotic
BD_p_jsu_gjrgarch= BD_jsu_gjrgarch$pvalue_twosided_asymptotic
BD_p_jsu_carr= BD_jsu_carr$pvalue_twosided_asymptotic

# FZL Function
FZL_sstd_nor_cop= mean(FZLoss(actual, VaR_sstd_nor_cop, ES_sstd_nor_vol, alpha=0.05))
FZL_sstd_t_cop= mean(FZLoss(actual, VaR_sstd_t_cop, ES_sstd_t_vol, alpha=0.05))
FZL_sstd_garch= mean(FZLoss(actual, VaR_sstd_garch, ES_sstd_garch, alpha = 0.05))
FZL_sstd_gjrgarch= mean(FZLoss(actual, VaR_sstd_gjrgarch, ES_sstd_gjrgarch,alpha = 0.05))
FZL_sstd_nagarch= mean(FZLoss(actual, VaR_sstd_nagarch, ES_sstd_nagarch,alpha = 0.05))
FZL_sstd_carr= mean(FZLoss(actual, VaR_sstd_carr, ES_sstd_carr, alpha=0.05))
FZL_msgarch= mean(FZLoss(actual, bid_msgarch_var_05, bid_msgarch_es_05, alpha=0.05))
FZL_gas= mean(FZLoss(actual, bid_gas_var_05, bid_gas_es_05, alpha=0.05))
FZL_robgarch= mean(FZLoss(actual, bid_rob_garch_var_05, bid_rob_garch_es_05, alpha=0.05))
FZL_ged_nor_cop= mean(FZLoss(actual, VaR_ged_nor_cop, ES_ged_nor_vol, alpha=0.05))
FZL_ged_t_cop= mean(FZLoss(actual, VaR_ged_t_cop, ES_ged_t_vol, alpha=0.05))
FZL_ged_garch= mean(FZLoss(actual, VaR_ged_garch, ES_ged_garch, alpha = 0.05))
FZL_ged_gjrgarch= mean(FZLoss(actual, VaR_ged_gjrgarch, ES_ged_gjrgarch,alpha = 0.05))
FZL_ged_nagarch= mean(FZLoss(actual, VaR_ged_nagarch, ES_ged_nagarch,alpha = 0.05))
FZL_ged_carr= mean(FZLoss(actual, VaR_ged_carr, ES_ged_carr, alpha=0.05))

FZL_sged_nor_cop= mean(FZLoss(actual, VaR_sged_nor_cop, ES_sged_nor_vol, alpha=0.05))
FZL_sged_t_cop= mean(FZLoss(actual, VaR_sged_t_cop, ES_sged_t_vol, alpha=0.05))
FZL_sged_garch= mean(FZLoss(actual, VaR_sged_garch, ES_sged_garch, alpha = 0.05))
FZL_sged_gjrgarch= mean(FZLoss(actual, VaR_sged_gjrgarch, ES_sged_gjrgarch,alpha = 0.05))
FZL_sged_nagarch= mean(FZLoss(actual, VaR_sged_nagarch, ES_sged_nagarch,alpha = 0.05))
FZL_sged_carr= mean(FZLoss(actual, VaR_sged_carr, ES_sged_carr, alpha=0.05))
FZL_jsu_nor_cop= mean(FZLoss(actual, VaR_jsu_nor_cop, ES_jsu_nor_vol, alpha=0.05))
FZL_jsu_t_cop= mean(FZLoss(actual, VaR_jsu_t_cop, ES_jsu_t_vol, alpha=0.05))
FZL_jsu_garch= mean(FZLoss(actual, VaR_jsu_garch, ES_jsu_garch, alpha = 0.05))
FZL_jsu_gjrgarch= mean(FZLoss(actual, VaR_jsu_gjrgarch, ES_jsu_gjrgarch,alpha = 0.05))
FZL_jsu_nagarch= mean(FZLoss(actual, VaR_jsu_nagarch, ES_jsu_nagarch,alpha = 0.05))
FZL_jsu_carr= mean(FZLoss(actual, VaR_jsu_carr, ES_jsu_carr, alpha=0.05))

AE = rbind(AE_ged_nor_cop,AE_ged_t_cop,AE_ged_garch,AE_ged_nagarch,AE_ged_gjrgarch,AE_ged_carr,AE_msgarch,AE_robgarch,AE_gas,AE_sstd_nor_cop,AE_sstd_t_cop,AE_sstd_garch,AE_sstd_nagarch,AE_sstd_gjrgarch,AE_sstd_carr,AE_sged_nor_cop,AE_sged_t_cop,AE_sged_garch,AE_sged_nagarch,AE_sged_gjrgarch,AE_sged_carr,AE_jsu_nor_cop,AE_jsu_t_cop,AE_jsu_garch,AE_jsu_nagarch,AE_jsu_gjrgarch,AE_jsu_carr)
UC= rbind(uc_ged_nor_cop,uc_ged_t_cop,uc_ged_garch,uc_ged_nagarch,uc_ged_gjrgarch,uc_ged_carr,uc_msgarch,uc_robgarch,uc_gas,uc_sstd_nor_cop,uc_sstd_t_cop,uc_sstd_garch,uc_sstd_nagarch,uc_sstd_gjrgarch,uc_sstd_carr,uc_sged_nor_cop,uc_sged_t_cop,uc_sged_garch,uc_sged_nagarch,uc_sged_gjrgarch,uc_sged_carr,uc_jsu_nor_cop,uc_jsu_t_cop,uc_jsu_garch,uc_jsu_nagarch,uc_jsu_gjrgarch,uc_jsu_carr)
CC= rbind(cc_ged_nor_cop,cc_ged_t_cop,cc_ged_garch,cc_ged_nagarch,cc_ged_gjrgarch,cc_ged_carr,cc_msgarch,cc_robgarch,cc_gas,cc_sstd_nor_cop,cc_sstd_t_cop,cc_sstd_garch,cc_sstd_nagarch,cc_sstd_gjrgarch,cc_sstd_carr,cc_sged_nor_cop,cc_sged_t_cop,cc_sged_garch,cc_sged_nagarch,cc_sged_gjrgarch,cc_sged_carr,cc_jsu_nor_cop,cc_jsu_t_cop,cc_jsu_garch,cc_jsu_nagarch,cc_jsu_gjrgarch,cc_jsu_carr)
DQ= rbind(DQ_ged_nor_cop,DQ_ged_t_cop,DQ_ged_garch,DQ_ged_nagarch,DQ_ged_gjrgarch,DQ_ged_carr,DQ_msgarch,DQ_robgarch,DQ_gas,DQ_sstd_nor_cop,DQ_sstd_t_cop,DQ_sstd_garch,DQ_sstd_nagarch,DQ_sstd_gjrgarch,DQ_sstd_carr,DQ_sged_nor_cop,DQ_sged_t_cop,DQ_sged_garch,DQ_sged_nagarch,DQ_sged_gjrgarch,DQ_sged_carr,DQ_jsu_nor_cop,DQ_jsu_t_cop,DQ_jsu_garch,DQ_jsu_nagarch,DQ_jsu_gjrgarch,DQ_jsu_carr)
McF= rbind(McF_p_ged_nor_cop,McF_p_ged_t_cop,McF_p_ged_garch,McF_p_ged_nagarch,McF_p_ged_gjrgarch,McF_p_ged_carr,McF_p_msgarch,McF_p_robgarch,McF_p_gas,McF_p_sstd_nor_cop,McF_p_sstd_t_cop,McF_p_sstd_garch,McF_p_sstd_nagarch,McF_p_sstd_gjrgarch,McF_p_sstd_carr,McF_p_sged_nor_cop,McF_p_sged_t_cop,McF_p_sged_garch,McF_p_sged_nagarch,McF_p_sged_gjrgarch,McF_p_sged_carr,McF_p_jsu_nor_cop,McF_p_jsu_t_cop,McF_p_jsu_garch,McF_p_jsu_nagarch,McF_p_jsu_gjrgarch,McF_p_jsu_carr)
NF= rbind(NF_p_ged_nor_cop,NF_p_ged_t_cop,NF_p_ged_garch,NF_p_ged_nagarch,NF_p_ged_gjrgarch,NF_p_ged_carr,NF_p_msgarch,NF_p_robgarch,NF_p_gas,NF_p_sstd_nor_cop,NF_p_sstd_t_cop,NF_p_sstd_garch,NF_p_sstd_nagarch,NF_p_sstd_gjrgarch,NF_p_sstd_carr,NF_p_sged_nor_cop,NF_p_sged_t_cop,NF_p_sged_garch,NF_p_sged_nagarch,NF_p_sged_gjrgarch,NF_p_sged_carr,NF_p_jsu_nor_cop,NF_p_jsu_t_cop,NF_p_jsu_garch,NF_p_jsu_nagarch,NF_p_jsu_gjrgarch,NF_p_jsu_carr)
BD= rbind(BD_p_ged_nor_cop,BD_p_ged_t_cop,BD_p_ged_garch,BD_p_ged_nagarch,BD_p_ged_gjrgarch,BD_p_ged_carr,BD_p_msgarch,BD_p_robgarch,BD_p_gas,BD_p_sstd_nor_cop,BD_p_sstd_t_cop,BD_p_sstd_garch,BD_p_sstd_nagarch,BD_p_sstd_gjrgarch,BD_p_sstd_carr,BD_p_sged_nor_cop,BD_p_sged_t_cop,BD_p_sged_garch,BD_p_sged_nagarch,BD_p_sged_gjrgarch,BD_p_sged_carr,BD_p_jsu_nor_cop,BD_p_jsu_t_cop,BD_p_jsu_garch,BD_p_jsu_nagarch,BD_p_jsu_gjrgarch,BD_p_jsu_carr)
QL= rbind(QL_ged_nor_cop,QL_ged_t_cop,QL_ged_garch,QL_ged_nagarch,QL_ged_gjrgarch,QL_ged_carr,QL_msgarch,QL_robgarch,QL_gas,QL_sstd_nor_cop,QL_sstd_t_cop,QL_sstd_garch,QL_sstd_nagarch,QL_sstd_gjrgarch,QL_sstd_carr,QL_sged_nor_cop,QL_sged_t_cop,QL_sged_garch,QL_sged_nagarch,QL_sged_gjrgarch,QL_sged_carr,QL_jsu_nor_cop,QL_jsu_t_cop,QL_jsu_garch,QL_jsu_nagarch,QL_jsu_gjrgarch,QL_jsu_carr)
FZL= rbind(FZL_ged_nor_cop,FZL_ged_t_cop,FZL_ged_garch,FZL_ged_nagarch,FZL_ged_gjrgarch,FZL_ged_carr,FZL_msgarch,FZL_robgarch,FZL_gas,FZL_sstd_nor_cop,FZL_sstd_t_cop,FZL_sstd_garch,FZL_sstd_nagarch,FZL_sstd_gjrgarch,FZL_sstd_carr,FZL_sged_nor_cop,FZL_sged_t_cop,FZL_sged_garch,FZL_sged_nagarch,FZL_sged_gjrgarch,FZL_sged_carr,FZL_jsu_nor_cop,FZL_jsu_t_cop,FZL_jsu_garch,FZL_jsu_nagarch,FZL_jsu_gjrgarch,FZL_jsu_carr)

Table= cbind(AE,UC,CC,DQ,McF,NF,BD,QL,FZL)

write.csv(Table, "Table_05.csv")

# Similarly estime VaR and ES at 1% and 2.5% risk levels
