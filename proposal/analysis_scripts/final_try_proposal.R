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
climate_window <- window(climate_ts,start=start(spruce_window),end=end(spruce_window))
plotc(spruce_window)

log_spruce_window <- log(spruce_window)
t <- seq(1,length(spruce_window))
t2 <- t^2

d <- data.frame(cbind(log_spruce_window,t,t2))

#log_t <- log(t)
#m1 <- lm(log_spruce_window ~ log_t)
#summary(m1)

m2 <- lm(log_spruce_window ~ t+t2,data=d)
summary(m2)

i <- (m2$residuals/m2$fitted.values)
plotc(i)
acf(i)
pacf(i)

kpss.test(i)
adf.test(i)
pp.test(i)

model_ar2 <- arima(i,order=c(2,0,0))
summary(model_ar2)

model_arma11 <- arima(i,order=c(1,0,1))
