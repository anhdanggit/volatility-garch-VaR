setwd("E:/TSE/Financial Econometrics")
df <- read.delim("E:/TSE/Financial Econometrics/daily_rets_day_ret_RV_BV_RJ_2.txt", 
                 header=FALSE)
colnames(df) <- c("Date", "Return", "Volatility", "V4", "V5")
df$Return2 <- df$Return^2
summary(df)

#### 1. Descriptive Statistics ####
library(moments)
library(ggplot2)
library(zoo)

# Distribution
skewness(df$Return) # Return
kurtosis(df$Return)
qplot(df$Return, geom = 'histogram', binwidth = 0.008) + xlab('Return')

skewness(df$Return2) # Squared Return
kurtosis(df$Return2)
qplot(df$Return2, geom = 'histogram', binwidth = 0.0003) + xlab('Squared Return')

skewness(df$Volatility) # Realized Volatility
kurtosis(df$Volatility)
qplot(df$Volatility, geom = 'histogram', binwidth = 0.0001) + xlab('Realized Volatility')

# Plotting (Sample)
acf(df$Return, plot = TRUE, main = "ACF of Return with several lags") # Return
acf(df$Return2, plot = TRUE, main = "ACF of Squared Return with several lags") # Squared Return
acf(df$Volatility, plot = TRUE, main = "ACF of RV with several lags") # Realized Volatility

plot(as.ts(df$Return), main = "Daily Return", col = 4, ylab = "Return") # Return
plot(as.ts(df$Return2), main = "Daily Squared Return", col = 4, ylab = "Squared Return")
plot(as.ts(df$Volatility), main = "Realized Volatility", col = 4, ylab = "RV")

# Plotting (Rolling)
roll <- zoo(df)
window <- 120

plot(na.omit(rollapply(roll$Return, window, function(x) mean(x))), main = "Rolling Mean", xlab = "Time", ylab = "", col = 4)
plot(na.omit(rollapply(roll$Return, window, function(x) var(x))), main = "Rolling Variance", xlab = "Time", ylab = "", col = 4)
plot(na.omit(rollapply(roll$Return, window, function(x) skewness(x))), main = "Rolling Skewness", xlab = "Time", ylab = "", col = 4)
plot(na.omit(rollapply(roll$Return, window, function(x) kurtosis(x))), main = "Rolling Kurtosis", xlab = "Time", ylab = "", col = 4)


plot(na.omit(rollapply(roll$Volatility, window, function(x) mean(x))), main = "Rolling Mean", xlab = "Time", ylab = "", col = 2)
plot(na.omit(rollapply(roll$Volatility, window, function(x) var(x))), main = "Rolling Variance", xlab = "Time", ylab = "", col = 2)
plot(na.omit(rollapply(roll$Volatility, window, function(x) skewness(x))), main = "Rolling Skewness", xlab = "Time", ylab = "", col = 2)
plot(na.omit(rollapply(roll$Volatility, window, function(x) kurtosis(x))), main = "Rolling Kurtosis", xlab = "Time", ylab = "", col = 2)


#### 2. Fitting GARCH model ####

# (8) Riskmetriks ####
library(rugarch)
riskmetrik.spec = ugarchspec(mean.model=list(armaOrder=c(1,1), include.mean=TRUE),
                             variance.model=list(model="iGARCH"),
                             distribution.model = "ged",
                             fixed.pars = list(omega = 0))


riskmetrik.model = ugarchfit(spec = riskmetrik.spec, data = spyder.daily, out.sample = 500)

show(riskmetrik.model)

for.riskmetrik = ugarchforecast(riskmetrik.model, data=NULL, n.head=1, n.roll = 499, out.sample = 500)


# Out-of-sample

for.risk.var = as.ts(for.riskmetrik@forecast$sigmaFor[1,])^2

for.egarch.var = as.ts(for.model_2@forecast$sigmaFor[1,])^2

for.rv = as.ts(window(dat$RV, start = 1998))

# In sample

risk.var = as.ts(sigma(riskmetrik.model))^2

egarch.var = as.ts(sigma(model_2))^2

# rv = as.ts(dat$RV[1:1997])
rv = window(as.ts(dat$RV), end = 1997)

# Comparision of sigma between GARCH, eGARCH and iGARCH

# Out-of-sample

# Time series
plot(for.egarch.var, col=4, lty=2, main = "Out-of-sample Variance Estimation", 
     ylab = "Est. Var")
lines(for.rv, col=3, lty=1)
lines(for.risk.var, col=2, lty=1)
l1 = as.expression("GARCH")
l2 = as.expression("RV")
l3 = as.expression("RiskMetrics")
legend("topright", c(l1, l2, l3), col=c(4,2,3), lty=c(2,1,1), bty="n", cex=0.7)


# Boxplot
boxplot(for.risk.var, for.egarch.var, for.rv, names = c("RiskMetrics", "GARCH", "RV"), 
        col = c(2, 4, 3), ylab = "Est. Var", 
        main = "Boxplot for Estimated Variance Comparison (Out-of-sample)")


##### AR estimation for RV ####
library(FitAR)
library(forecast)

# Model selection
selection <- SelectModel(as.ts(df$Volatility[1:1997]), lag.max = 14, ARModel = "AR",
                         Criterion = "AIC", Best = 6); selection

# Fitting AR(12)
fit1 <- arima(rv, order = c(12,0,0))
fit2 <- Arima(as.ts(dat$RV), model = fit1)
fc <- window(fitted(fit2), start = 1998)

# Testing
print(Box.test(fit2$residuals, lag=1, type = "Ljung-Box"))
print(Box.test(fit2$residuals, lag=5, type = "Ljung-Box"))
print(Box.test(fit2$residuals, lag=9, type = "Ljung-Box"))

print(Box.test(fit2$residuals^2, lag=1, type = "Ljung-Box"))
print(Box.test(fit2$residuals^2, lag=8, type = "Ljung-Box"))
print(Box.test(fit2$residuals^2, lag=14, type = "Ljung-Box"))

print(archTest(fit2$residuals, lag = 4))
print(archTest(fit2$residuals, lag = 6))
print(archTest(fit2$residuals, lag = 8))


# Calculate RMSE and MAE
error_ar <- (dat$RV[1998:length(dat$RV)] - fc)
n_ar = length(dat$RV[1998:length(dat$RV)])
rmse = sqrt(sum(error_ar^2) / n_ar)
mae = sum(abs(error_ar)) / n_ar

# Plotting
plot(for.rv, col = 4, main = "Out-of-sample forecast with AR(12)",
     ylab = "var", lty = 1)
lines(fc, col = 2, lty = 1)
legend("topright", c("Actual", "Prediction"), col=c(4,2), lty=c(1,2), bty="n", cex=0.7)