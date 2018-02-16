setwd("E:/TSE/Financial Econometrics")
df <- read.delim("E:/TSE/Financial Econometrics/daily_rets_day_ret_RV_BV_RJ_2.txt", 
                 header=FALSE)
colnames(df) <- c("Date", "Return", "Volatility", "V4", "V5")
df$Return2 <- df$Return^2
summary(df)

#### 1. Descriptive Statistics ####
library(moments)
library(ggplot2)

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

# ACF and PACF
acf(df$Return, plot = TRUE) # Return
pacf(df$Return, plot = TRUE)

acf(df$Return2, plot = TRUE) # Squared Return
pacf(df$Return2, plot = TRUE)

acf(df$Volatility, plot = TRUE) # Realized Volatility
pacf(df$Volatility, plot = TRUE)


#### 2. Fitting GARCH model ####
library(fGarch)
gmodel <- garchFit(formula = ~arma(1,1) + garch(2,1), 
                   data = df$Return, trace = FALSE )
summary(gmodel)
plot(gmodel)
