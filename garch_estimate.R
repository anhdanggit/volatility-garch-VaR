
library(rugarch)
library(dplyr)
library(car)
library(timeSeries)
library(tseries)
library(lmtest)
library(xts)
library(stargazer)
library(MTS)

# (1) Load data ####

dat = read.table("data/daily_rets_day_ret_RV_BV_RJ_2.txt", row.names = 1,
                  col.names = c("t", "daily", "RV", "BV", "RJ"))
dat = dat[order(row.names(dat)), ]

spyder.daily = as.ts(dat$daily)
plot(spyder.daily)


# (2) GARCH Modeling ####

# order of model  -> choose by the AIC

## parameters:

### garch.model = "GARCH", "TGARCH", "apARCH", "eGARCH"
### garch.order = c(1,1), c(2,1), etc.
### dist = "ged", "norm", "std"
### y = series to estimate
### optional: arma.order, out.sample

model_garch = function (garch.model, garch.order, dist, 
                        arma.order = c(1,1), y = spyder.daily, out.sample = 500) {
  
  if (dist %in% c("std", "norm", "ged") != TRUE) {
    stop("Your input distribution model is not supported. Try: 'std', 'norm', 'ged'")
  }
  
  if (garch.model == "GARCH" | garch.model == "TGARCH") {
    # model specification
    model.spec = ugarchspec(variance.model = list(model = "fGARCH", garchOrder=garch.order, submodel = garch.model),
                              mean.model = list(include.mean=TRUE, armaOrder=arma.order), #fixed.pars = list(alpha1=0),
                              distribution.model = dist) # "norm", "std", "ged": for student
    # model fitting
    fit.model = ugarchfit(spec = model.spec, data = y, out.sample = out.sample) # out-of-sample
  } 
  
  else if (garch.model == "eGARCH" | garch.model == "apARCH") {
    model.spec = ugarchspec(variance.model = list(model = garch.model, garchOrder=garch.order),
                            mean.model = list(include.mean=TRUE, armaOrder=arma.order), 
                            distribution.model = dist) 
    # model fitting
    fit.model = ugarchfit(spec = model.spec, data = y, out.sample = out.sample) 
  }
  else {
    stop("Your input GARCH model is not supported. Try: GARCH, TGARCH, eGARCH, apARCH ")
  }
  
  return(fit.model)
}

# (3) GARCH Forecast ####

forecast_garch = function(fit.model, data=NULL, n.head=1, n.roll=499, out.sample=500) {
  for.model = ugarchforecast(fit.model, data=NULL, n.head=1, n.roll = 499, out.sample = 500)
}


# (5) out-of-sample Diagnostic ####

accuracy.vol =  function(for.model, 
                     realize = dat$RV[1998:length(dat$RV)]){
  n = length(realize)
  
  sigma.hat = for.model@forecast$sigmaFor[1,]
  error = realize - sigma.hat^2 
  
  RMSE = sqrt(sum(error^2) / n)
  MAE = sum(abs(error)) / n
  
  sol = list(RMSE = RMSE, MAE = MAE)
  return(sol)
}

accuracy.return = function(for.model,
                           realize = dat$daily[1998:length(dat$daily)]) {
  n = length(realize)
  error = (dat$daily[1998:length(dat$RV)] - for.model@forecast$seriesFor[1,]) 
  
  RMSE = sqrt(sum(error^2) / n) 
  MAE = sum(abs(error)) / n
  
  sol = list(RMSE = RMSE, MAE = MAE)
  return(sol)
}

## Compute std. residuals

std.res = function (for.model, y = spyder.daily[1998:length(spyder.daily)]) {
  
  yhat = for.model@forecast$seriesFor[1,]
  epsilon = y - yhat
  sigma.hat = for.model@forecast$sigmaFor[1,]
  
  z = epsilon /sigma.hat
  
  return(z)
}

## Test residual
boxtest.z = function (for.model) {
  z = std.res(for.model)
  print("lag 1")
  print(Box.test(z, lag=1, type = "Ljung-Box"))
  print("lag 5")
  print(Box.test(z, lag=5, type = "Ljung-Box"))
  print("lag 9")
  print(Box.test(z, lag=9, type = "Ljung-Box"))
}


## Test sq-residual
boxtest.sqz = function (for.model) {
  z = std.res(for.model)
  print("lag 1")
  print(Box.test(z^2, lag=1, type = "Ljung-Box"))
  print("lag 8")
  print(Box.test(z^2, lag=8, type = "Ljung-Box"))
  print("lag 14")
  print(Box.test(z^2, lag=14, type = "Ljung-Box"))
}

## Test ARCH
archtest.z = function (for.model) {
  z = std.res(for.model)
  print("lag 4")
  print(archTest(z, lag=4))
  print("lag 6")
  print(archTest(z, lag=6))
  print("lag 8")
  print(archTest(z, lag=8))
}




# (3) VaR #####

## (3.1) Hit Sequence ####


hit.seq.insample = function(fit.model, alpha = 0.05) {
  VaR = quantile(fit.model, alpha)
  # Hit Sequence
  n = length(spyder.daily[1:1997])
  hit = rep(0, n)
  for (i in 1:1997){
    if(spyder.daily[i] < VaR[i]){
      hit[i] = 1
    }
  }
  
  return(hit)
}

VaR.outsample = function(for.model, dist, alpha = 0.05, 
                             mu=0, sigma=1, lambda=-0.5, skew=1, shape=5) {
  critical = qdist(distribution = dist , p = alpha, mu=mu, sigma=sigma, lambda = lambda, skew=skew, shape=shape)
  
  VaR = for.model@forecast$seriesFor[1,] - abs(critical*for.model@forecast$sigmaFor[1,])
  return(VaR)
  }


hit.seq.outsample = function(for.model, dist, alpha = 0.05) {
  VaR = VaR.outsample(for.model, dist, alpha)
  # Hit Sequence
  n = length(spyder.daily[1998:length(spyder.daily)])
  hit = rep(0, n)
  for (i in 1:n){
    if(spyder.daily[i] < VaR[i]){
      hit[i] = 1
    }
  }
  return(hit)
}

violation.diff = function(hit, fit.model, alpha = 0.05){
  VaR = quantile(fit.model, alpha)
  n = length(hit)
  diff=rep(0,n)
  diff[which(hit==1)] = spyder.daily[which(hit==1)] - VaR[which(hit==1)]
  return(diff)
}



## (3.2) Test1: Uncondtional Ber(alpha) ####

test1.uncond = function(hit, alpha = 0.05){
  ### pi = T1 / T (observed)
  n = length(hit)
  T1 = length(which(hit==1))
  T0 = n - T1
  pi = T1 / n
  
  # likelihood of unrestricted (pi)
  L.pi = (1-pi)^T0 * pi^T1 
  
  ### alpha (restricted)
  L.alpha = (1-alpha)^T0 * alpha^T1
  
  ### Log-likelihood ratio test
  LR.unc = 2*(log(L.pi) - log(L.alpha))
  LR.unc
}

## (3.3) Test 2: Independence ####

test2.ind = function(hit, alpha = 0.05){
  # T00 
  n = length(hit)
  T00 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 0 & hit[i+1] == 0){
      T00 = T00 + 1
    }
  }
  # T01 
  T01 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 0 & hit[i+1] == 1){
      T01 = T01 + 1
    }
  }
  # T10 
  T10 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 1 & hit[i+1] == 0){
      T10 = T10 + 1
    }
  }
  # T11 
  T11 = 0
  for (i in 1: (n-1)){
    if (hit[i] == 1 & hit[i+1] == 1){
      T11 = T11 + 1
    }
  }
  
  # Likelihood of pi
  pi.01 = T01 / (T00 + T01)
  pi.11 = T11 / (T10 + T11)
  pi.00 = 1 - pi.01
  pi.10 = 1 - pi.11
  
  L.pi2 = (1 - pi.01)^T00 * pi.01^T01 * (1-pi.11)^T10 * pi.11^T11
  L.alpha2 = (1 - alpha)^T00 * alpha^T01 * (1-alpha)^T10 * alpha^T11
  LRcc = 2 * ( log(L.pi2) - log(L.alpha2) )
  LRcc
}


## (3.3) Test 2: Independence ####

ttest <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1]-val)/co[coefnum,2]
  pval <- 2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
  return(c(tstat, pval))
}

test3.explain = function(hit, set) {
  n = length(hit)
  T = length(dat$BV)-1
  if (set == "in-sample") {
    test3 = lm(hit[2:length(hit)] ~ dat$BV[1:1996] + dat$RJ[1:1996])
  } else {
    test3 = lm(hit[2:length(hit)] ~ dat$BV[1998:T] + dat$RJ[1998:T])
  }
  print("b0=alpha")
  print(ttest(test3, 1, 0.05))
  print("b1 = 0")
  print(ttest(test3, 2, 0))
  print("b2 = 0")
  print(ttest(test3, 3, 0))

}
  

## (7) Results ####

# (7.1) GARCH Estimate ####

# Generalized error
model_1 = model_garch(garch.model = "GARCH", garch.order = c(2,1), dist = "ged")
model_2 = model_garch(garch.model = "eGARCH", garch.order = c(2,1), dist = "ged")
model_2B = model_garch(garch.model = "TGARCH", garch.order = c(1,1), dist = "ged")
model_2C = model_garch(garch.model = "apARCH", garch.order = c(1,1), dist = "ged")

# norm error
model_3 = model_garch(garch.model = "GARCH", garch.order = c(2,1), dist = "norm")
model_4 = model_garch(garch.model = "eGARCH", garch.order = c(2,1), dist = "norm")

# student error
model_5 = model_garch(garch.model = "GARCH", garch.order = c(2,1), dist = "std")
model_6 = model_garch(garch.model = "eGARCH", garch.order = c(2,1), dist = "std") 


# (7.1) GARCH Forecast ####
for.model_1 = forecast_garch(model_1)
for.model_2 = forecast_garch(model_2)

for.model_3 = forecast_garch(model_3)
for.model_4 = forecast_garch(model_4)

for.model_5 = forecast_garch(model_5)
for.model_6 = forecast_garch(model_6)


# (7.2) Out-sample Diagnostics ####

# Apply model!!!!!-------#
model = model_1
for.model = for.model_1
#------------------------#

accuracy.vol(for.model)
boxtest.z(for.model)
boxtest.sqz(for.model)
archtest.z(for.model)
jarque.bera.test(std.res(for.model))


# (7.3) Plot Models ####

# Apply model!!!!!-------#
model = model_1
for.model = for.model_1
#------------------------#

show(model)
plot(model)

plot(for.model)

# (7) VaR Test ####

# Apply model!!!!!-------#
model = model_1
for.model = for.model_1
#------------------------#

hit = hit.seq.insample(model)
length(which(hit==1)) / length(hit) # violation rate
diff = violation.diff(hit, model) # violation diff.

#plot.ts(diff, ylim = c(-0.04, 0))
test1.uncond(hit)
test2.ind(hit)
test3.explain(hit, set = "in-sample")


hit.out = hit.seq.outsample(for.model, dist = "ged")
length(which(hit1.out==1)) / length(hit.out) # violation rate
diff.out = violation.diff(hit, model) # violation diff.

#plot.ts(diff.out)
test1.uncond(hit.out)
test2.ind(hit.out)
test3.explain(hit.out, set = "out-sample")

