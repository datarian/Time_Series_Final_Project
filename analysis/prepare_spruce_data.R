library(itsmr)
library(tseries)

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

spruce_window <- window(spruce_sup_900_ts,start=1400,end=1800)

plot(spruce_sup_900_ts,
     main="Spruce, > 900 m",
     ylab="tree ring width [1/100 mm]",
     xlab="year")

rwl <- readRDS("data/rwl_900+.Rds")

rwl_mean <- apply(rwl,1,mean,na.rm=T)
rwl_mean_ts <- ts(rwl_mean,end=2017)
rwl_window <- window(rwl_mean_ts,start=1400,end=1800)

plot(spruce_window)
