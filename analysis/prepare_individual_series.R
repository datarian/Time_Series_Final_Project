#### Individual series
# This script parses the raw data file 900+.txt and builds a dataframe in 'rwl'
# format, a standard for dendrochronological data.
# It also computes the mean value chronology and does som preliminary analysis
# of the data.

library(readtext)
library(stringr)
library(dplR)
library(itsmr)
library(MASS)

file_contents <- readtext("data/900+.txt")

# Extract the individual data series from the unstructured text file
pattern <- "\\sNr\\. ([0-9]*)\\.0 Dat: ([0-9]{3,4}).*\\n.*\\n.*\\n{1,2}((?:[ ]{1,3}[0-9 -]*\\n){3,})"
parsed <- str_match_all(file_contents,pattern)
rm(file_contents) # Free up memory

# Extract the sample numbers and their corresponding dated year
chronology <- as.data.frame(parsed[[1]][,-1],stringsAsFactors = F)
rm(parsed) # Free up memory
chron_numbers <- as.numeric(chronology[,1])
chron_years <- as.numeric(chronology[,2])

# Extract the ring width sequences for individual samples.
pat <- "\\s{1,3}[0-9 -]{8,11}\\s+((([0-9]{1,})\\s*)+)\\\n"
rings <- str_match_all(chronology[,3],pat)

# We end up with a list of matrices containing the group matches from regex.
# We need the 2nd column of each matrix which contains the data. Extract and
# coerce to a vector:
chron_rings <- list()
for (i in 1:nrow(chronology)) {
    chron_rings[[i]] <- scan(text = rings[[i]][,2])
}

# We now have three objects to construct the aggregated chronology dataframe:
# chron_numbers: the sample numbers (identifier)
# chron_years: the dated years of the samples
# chron_rings: The sequence of yearly ring widths of the samples

# construct the sequence of years spanning the data. For the oldest year, we need
# the minimum of dated year - age of tree
oldest <- 10000 # initialize with a ridiculously high number
for(i in 1:length(chron_years)){
    currentsample <- chron_years[i]-(length(chron_rings[[i]])-1)
    oldest <- ifelse(currentsample < oldest, currentsample, oldest)
}

newest <- max(chron_years)
year <- seq(oldest,newest,by=1)

# Result dataframe. Structure: First column is the year, subsequent columns
# represent one sample each, colname is the sample number. Years without
# tree ring widhts are coded as NA
rwl_df <- as.data.frame(year); rownames(rwl_df) <- year

for(i in 1:length(chron_numbers)){
    datedyear <- chron_years[i]
    fillbefore <- numeric()
    fillafter <- numeric()
    if(datedyear > oldest){
        fillbefore <- rep(NA,times=(datedyear-(length(chron_rings[[i]])-1)-oldest))
    }
    if(datedyear<newest){
        fillafter <- rep(NA,times=newest-datedyear)
    }
    rings <- chron_rings[[i]]
    series <- as.data.frame(c(fillbefore,rings,fillafter))
    colnames(series) <- chron_numbers[i]
    rwl_df <- cbind(rwl_df,series)
}

# remove year column:
rwl_df <- rwl_df[,-1]

# Some prior analysis:
dplR::rwl.report(rwl_df)
dplR::summary.rwl(rwl_df)


rwi <- dplR::detrend(rwl_df,make.plot = F,method = "Mean")
# plot the obtained chronology. A nice feature is the information on the sample
# depth!
rwi.crn <- dplR::chron(rwi)
plot(rwi.crn,add.spline=T,nyrs=20)



#*******************************************************************************
# Manually construct the time series from the yearly mean ring widhts, then
# dividing by the SD to remove some homoscedasticity
mean_rwl <- apply(rwl_df,1,mean,na.rm=T)
sd_rwl <- apply(rwl_df,1,sd,na.rm=T)

depth <- function(x){
    sum(!is.na(x))
}

depth_rwl <- apply(rwl_df,1,depth)


# Creating initial ts
rwl_ts <- ts(mean_rwl,end=2017)
rwl_ts_stab <- rwl_ts/sd_rwl
plotc(rwl_ts)


# limiting the series to the window 1400 ... 1800 (start chosen because of very low
# sample depth before that year)
rwl_ts_window <- window(rwl_ts_stab,start=c(1400),end=c(1800))
t <- seq(1,length(rwl_ts_window))
plotc(rwl_ts_window)



#*************************************
# First approach: Box-Cox transform and differencing at lag 1
# time


# Box-Cox-Transform:
lambdas <- boxcox(rwl_ts_window~t)
l <- lambdas$x[which.max(lambdas$y)] # this is the MLE lambda to transform data

rwl_ts_window_boxcox <- (rwl_ts_window^l-1)/l # Box-Cox transformation

# differencing:
rwl_ts_window_boxcox_diff <- diff(rwl_ts_window_boxcox)
plotc(rwl_ts_window_boxcox_diff)


acf(rwl_ts_window_boxcox_diff)
pacf(rwl_ts_window_boxcox_diff)



#******************************************
# Second approach: Log transform, polynomials
rwl_ts_window_log <- log(rwl_ts_window)
plotc(rwl_ts_window_log)

t <- seq(1,length(rwl_ts_window))
t2 <- t^2

m <- lm(rwl_ts_window_log ~ t)
summary(m)

m2 <- lm(rwl_ts_window_log ~ t+t2)
summary(m2)
# Nope, order 1 seems okay.


acf(m$residuals) # looks exponential now.
pacf(m$residuals)


# Maybe another approach? -> Weighing by sample depth?
weight_function <- function(x){
    n <- sum(!is.na(x))
    1 - ((ncol(rwl_df)-n)/ncol(rwl_df))
}
weight_rwl <- apply(rwl_df, 1, weight_function)
