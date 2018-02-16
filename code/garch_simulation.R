library(timeSeries)
library(ggplot2)
library(tidyr)
library(psych)

## (0) Setting simulation parameters ###############

R = 100 # num of replicate
n = 1000 # num of obs

#(mu, omega, alpha, beta, theta)
design = 1

#equation of ht
eq = 1 

#distribution of zt
dist = "norm" 

# distribution used in MLE
garch.dist = dist

#---------------------------------------------------#
# The resuls is in part 3 
if (design == 1) {
  mu = 0 
  omega = 0.01
  alpha = 0.05
  beta = 0.9
  theta = 0
}

if (design == 2) {
  mu = 0 
  omega = 0.01
  alpha = 0.05
  beta = 0.6
  theta = 2
}

if (dist == "norm"){
  z.random = function(n){
    rnorm(n,0,1)
  }
}

if (dist == "std"){
  z.random = function(n){
    rt(n, df=10)
  }
}

if (eq == 1){
  h.equation = function(omega, alpha, Mean, z, series, theta){
    omega + alpha*c(Mean, (z[-length(series)] - theta)^2)
  }
}
if (eq == 2){
  h.equation = function(omega, alpha, Mean, z, series, theta){
    omega + alpha*c(Mean, (z[-length(series)] - theta*sqrt(Mean))^2)
  }
}

  
#---------------------------------------#
## (2) Create the GARCH series #######

garch.sim <- function (n, n.start = 250, mu, omega, alpha, beta, theta, distribution=dist) {
  
  n = n + n.start # burnt-in
  mu = mu; omega = omega; alpha = alpha; beta = beta; theta = theta
  
  z = z.random(n)
  
  epsi = rep(0, n)
  h = rep(0, n)
  r = rep(0, n)
  
  for (i in 2:n) {
    h[i] = omega + alpha * (epsi[i-1] + mu - theta)^2 + beta * h[i-1]
    epsi[i] = z[i] * sqrt(h[i])
    r[i] = mu + epsi[i]
  }
  r = r[(n.start+1):n]
  
  return(r)
  
}


#### (2) Step-by-step: GARCH by MLE ####


### (2.1) Setting initial values ####

garchInit = function(series) {
  Mean = mean(series); Var = var(series); S = 1e-6
  params = c(mu = Mean, omega = 0.1*Var, alpha = 0.1, beta = 0.8, theta = 0)
  lowerBounds = c(mu = -10*abs(Mean), omega = S^2, alpha = S, beta = S, theta = -3)
  upperBounds = c(mu = 10*abs(Mean), omega = 100*Var, alpha = 1-S, beta = 1-S, theta = 3)
  cbind(params=params, lowerBounds=lowerBounds, upperBounds=upperBounds)
}

### (2.2) z distribution ####

if (garch.dist == "norm"){
  garchDist = function(z, hh){dnorm(x = z/hh)/hh}
}

if (garch.dist == "std"){
  garchDist = function(z,hh){dt(x=z/hh, df=10)/hh}
}

### (2.3) LLH Obj function ####

garchLLH = function(parm, series){
  mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]; theta = parm[5]
  z = (series - mu) #
  Mean = mean((z^2)) # first-term
  
  # Use Filter Representation
  e = h.equation(omega, alpha, Mean, z, series, theta)
  
  h = filter(e, beta, "r", init = Mean)
  hh = sqrt(abs(h))
  llh = -sum(log(garchDist(z, hh))) # use garchDist
  llh
}

### (2.4) MLE Optimize ####
garch11Fit = function(series, params, lowerBounds, upperBounds){
  fit = nlminb(start = params, objective = garchLLH, #use garchLLH 
               lower = lowerBounds, upper = upperBounds,
               control = list(trace=10), series = series)
  fit
}


## (3) Results ###############

mat = matrix(0, nrow = R, ncol = 5)

for (i in 1:R){
  
  # create r1
  r = garch.sim(n = n, n.start = 250, 
                mu = mu, omega = omega, alpha = alpha, beta = beta, theta = theta) # change parameter
  r = ts(r)
  
  garch.init = garchInit(r)
  params = garch.init[,1]
  lowerBounds = garch.init[,2]
  upperBounds = garch.init[,3]
  
  fit = garch11Fit(r, params, lowerBounds, upperBounds)
  
  mat[i,] = fit$par
}

# estimate averaring over replicates
est = colMeans(mat) 
sd.est = colSds(mat) 
# reported results
est
sd.est

# visualize
if (n == 100){
  est.param.100 = data.frame(n = rep("100", R), mat)
  names(est.param.100) = c("n","mu", "omega", "alpha", "beta", "theta")
  est.param.100 = gather(est.param.100, parameter, estimation, mu:theta, factor_key=TRUE)
}

if (n==1000){
  est.param.1000 = data.frame(n = rep("1000", R), mat)
  names(est.param.1000) = c("n","mu", "omega", "alpha", "beta", "theta")
  est.param.1000 = gather(est.param.1000, parameter, estimation, mu:theta, factor_key=TRUE)
}
est.param.1000 = data.frame(n = rep("1000", R), mat)
names(est.param.1000) = c("n","mu", "omega", "alpha", "beta", "theta")

estimate.result = rbind(est.param.100, est.param.1000)

est.param = est.param.1000
ggplot(data=est.param.100, aes(x = mu)) +
  geom_histogram(bins = 30) +
  geom_vline(data = est.param.100, aes(xintercept = mean(mu)), linetype = "dashed", size = 1) + 
  geom_vline(aes(xintercept=0), color = "red", size = 1) 

ggplot(data=est.param.100, aes(x = omega)) +
  geom_histogram(bins = 30) +
  geom_vline(data = est.param.100, aes(xintercept = mean(omega)), linetype = "dashed", size = 1) + 
  geom_vline(aes(xintercept=0.01), color = "red", size = 1)

ggplot(data=est.param.100, aes(x = alpha)) +
  geom_histogram(bins = 30) +
  geom_vline(data = est.param.100, aes(xintercept = mean(alpha)), linetype = "dashed", size = 1) + 
  geom_vline(aes(xintercept=0.05), color = "red", size = 1) 

ggplot(data=est.param.100, aes(x = beta)) +
  geom_histogram(bins = 30) +
  geom_vline(data = est.param.100, aes(xintercept = mean(beta)), linetype = "dashed", size = 1) + 
  geom_vline(aes(xintercept=0.9), color = "red", size = 1)

ggplot(data=est.param.100, aes(x = theta)) +
  geom_histogram(bins = 30) +
  geom_vline(data = est.param.100, aes(xintercept = mean(theta)), linetype = "dashed", size = 1) + 
  geom_vline(aes(xintercept=0), color = "red", size = 1) 
