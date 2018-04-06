# Time series project: Dendro data 

library(treeclim)
library(dplR) # package containing the Arizona dendro dataset
library(tseries)

# Literature
# [1] https://cran.r-project.org/web/packages/dplR/dplR.pdf
# [2] http://waldwachstum.wzw.tum.de/fileadmin/publications/2012_Pretzsch_dvffa.pdf
# [3] https://www.rdocumentation.org/packages/treeclim/versions/2.0.0
# https://www.ncdc.noaa.gov/cag/statewide/time-series/2/pcp/1/2/1895-2018?base_prd=true&firstbaseyear=1901&lastbaseyear=2000 (US weather data: Temp.,
# precipitation, drought index and many more)
# [4] https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1583&context=jmasm (Variance stabilizing power transformation for time series)


# data [1], p.37: This data set includes ring-width measurements for ponderosa pine (Pinus ponderosa) increment cores collected at the Gus Pearson Natural Area (GPNA) 
# in northern Arizona, USA. There are 58 series from 29 trees (2 cores per tree). Data are further described by Biondi and Qeadan (2008) and references therein.
data(gp.rwl)

str(gp.rwl)
head(gp.rwl)
data <- as.matrix(gp.rwl)
means <- rowMeans(data, na.rm=T)

plot(y=means, x=1:length(means), type="l")

# differenced time sereis
data_diff <- diff(means)

plot(y=data_diff, x=1:length(data_diff), type="l")

# H0: stationary cannot be rejected. But data 
kpss.test(data_diff)

# is this MA(2)?
acf(data_diff)
pacf(data_diff)

# Artificial MA(2) process, quite similar to dendro process
set.seed(999)
noise <- rnorm(400, 1, 1)
n <- length(noise)
x <- noise[3:n] - 0.5*noise[2:(n-1)] - 0.2*noise[1:(n-2)]
plot(x, type="l")
acf(x)
pacf(x)

# -------------------------------------------------------------------

# Nice :))))))
data(anos1)
anos1
str(anos1)

yrs <- as.numeric(rownames(anos1))
# Plot rings with data of two radii from same individual tree
res <- plotRings(yrs, anos1[,4], trwW = anos1[,5],
                 species.name = "Cedrela odorata")


