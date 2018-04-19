library(itsmr)
library(tseries)
library(MASS)
library(forecast)

spruce_sup900 <- read.table("./data/62167_picea_abv900.txt",
                            skip=4,nrows=68)
spruce_sup900_lastrow <- read.table("./data/62167_picea_abv900.txt",
                                    skip=72,nrows=1)
ssup9 <- spruce_sup900[,-c(1,2,3)]
ssup9l <- spruce_sup900_lastrow[,-c(1,2,3)]

series_spruce_sup900 <- c()
for (i in 1:nrow(ssup9)) {
    for (j in 1:ncol(ssup9)) {
        series_spruce_sup900 <- c(series_spruce_sup900,as.integer(ssup9[i,j]))
    }
}
series_spruce_sup900 <- c(series_spruce_sup900,as.integer(ssup9l))
spruce_sup_900_ts <- ts(series_spruce_sup900, start=1332)


plot(spruce_sup_900_ts,
     main="Spruce, > 900 m",
     ylab="tree ring width [1/100 mm]",
     xlab="year")

climate <- read.csv("./data/loehle_global_temp_reconstruction.csv",header=T)
climate_ts <- ts(climate$Temp..Anom.,start=16,end=1935)


# getting windows of the two series from 1332 to 1800:
spruce_window <- window(spruce_sup_900_ts,end=c(1800))

log_spruce_window <- log(spruce_window)
t <- seq(1,length(spruce_window))
log_t <- log(t)

d <- data.frame(cbind(log_spruce_window,log_t))

t_spruce_window <- nls(spruce_window ~ alpha * t^beta * exp(t),
                       data = d,
                       start=list(alpha=0.4, beta=0.8)
                       )


model1 <- lm(log_spruce_window ~ log_t)

plot(resid(model1))

kpss.test(resid(model1))

acf(resid(model1))

pacf(resid(model1))

climate_window <- window(climate_ts,start=start(spruce_window),end=end(spruce_window))

plotc(log(spruce_window))

# First try: yearly temperature anomalies as a covariate.
t <- seq(1,length(spruce_window))
climate_model <- lm(spruce_window~t+climate_window)
summary(climate_model)

ar1_residuals <- diff(climate_model$residuals)
plotc(ar1_residuals)

acf(ar1_residuals)
pacf(ar1_residuals)

ar2_residuals <- diff(diff(climate_model$residuals))
plotc(ar2_residuals)

acf(ar2_residuals)
pacf(ar2_residuals)

ar(climate_model$residuals,method="mle")

Box.test(climate_model$residuals,lag=5,type="Ljung-Box")
# H0: The residuals are iid distributed
# H1: The residuals are not iid distributed

# --> No luck

###############

# Find lambdas. MLE-estimation
lambdas <- boxcox(spruce_window~t+offset_climate)
l <- lambdas$x[which.max(lambdas$y)] # this is the MLE lambda to transform data

spruce_window_boxcox <- (spruce_window^l-1)/l # Box-Cox transformation

plotc(spruce_window_boxcox)

# New approach: offset between climatic data and treerings?
climate_window_long <- window(climate_ts,start=(start(spruce_window)[1]-50),end=end(spruce_window))

offset <- numeric(50)
pvalue <- numeric(50)

# Fit different models for a negative offset of 50 ... 1 years
for(i in 50:1) {
    model <- lm(spruce_window_boxcox~t+window(climate_ts,start=(start(spruce_window)[1]-i+1),end=(end(spruce_window)[1]-i+1)))
    offset[(51-i)] <- i
    pvalue[(51-i)] <- summary(model)$coefficients[,4][3] # get the p-value for climate
}

plot(pvalue~offset) # the smallest p-value results for an offset of 49 years.
offset[which.min(pvalue)] # Starting from a 10-year offset, the climatic covariate becomes significant.
# 10 years are chosen because it is intuitively more likely that this influences tree growth compared to a longer offset

# construct new climate window + model
offset_climate <- window(climate_ts,start=(start(spruce_window)[1]-10),end=(end(spruce_window)[1]-10))
climate_offset_model <- lm(spruce_window~t+offset_climate)
summary(climate_offset_model)
plot(climate_offset_model$residuals)

# Finally, transform treering data first:
# Box-Cox transformation of the data to obtain normally-distributed rv:

# Find lambdas. MLE-estimation
lambdas <- boxcox(spruce_window~t+offset_climate)
l <- lambdas$x[which.max(lambdas$y)] # this is the MLE lambda to transform data

spruce_window_boxcox <- (spruce_window^l-1)/l # Box-Cox transformation

# New model combining climatic offset and transformed data
climate_offset_model_boxcox <- lm(spruce_window_boxcox~t+offset_climate)
summary(climate_offset_model_boxcox)
plotc(climate_offset_model_boxcox$residuals)
qqnorm(climate_offset_model_boxcox$residuals);qqline(climate_offset_model_boxcox$residuals)

plotc(diff(climate_offset_model_boxcox$residuals))
kpss.test(diff(climate_offset_model_boxcox$residuals))
Box.test(diff(climate_offset_model_boxcox$residuals),type="Ljung-Box") # independency okay
acf(diff(climate_model_boxcox$residuals))
pacf(diff(climate_model_boxcox$residuals))


## Trying the auto.arima method

if (require("forecast", character.only = TRUE)) {
    sup900.arima <- auto.arima(spruce_sup_900_ts, ic="bic")
    summary(sup900.arima)
    head(residuals(sup900.arima))
    coef(sup900.arima)
    acf(residuals(sup900.arima),plot=FALSE)
}


# Frequency analysis
redf.dat <- redfit(spruce_sup_900_ts, nsim = 1000)
par(tcl = 0.5, mar = rep(2.2, 4), mgp = c(1.1, 0.1, 0))
plot(redf.dat[["freq"]], redf.dat[["gxxc"]],
       ylim = range(redf.dat[["ci99"]], redf.dat[["gxxc"]]),
       type = "n", ylab = "Spectrum (dB)", xlab = "Frequency (1/yr)",
       axes = FALSE)
grid()
lines(redf.dat[["freq"]], redf.dat[["gxxc"]], col = "#1B9E77")
lines(redf.dat[["freq"]], redf.dat[["ci99"]], col = "#D95F02")
lines(redf.dat[["freq"]], redf.dat[["ci95"]], col = "#7570B3")
lines(redf.dat[["freq"]], redf.dat[["ci90"]], col = "#E7298A")
freqs <- pretty(redf.dat[["freq"]])
pers <- round(1 / freqs, 2)
axis(1, at = freqs, labels = TRUE)
axis(3, at = freqs, labels = pers)
mtext(text = "Period (yr)", side = 3, line = 1.1)
axis(2); axis(4)
legend("topright", c("dat", "CI99", "CI95", "CI90"), lwd = 2,
         col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
         bg = "white")
box()
par(op)
