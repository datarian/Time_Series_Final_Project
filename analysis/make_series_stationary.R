library(ggplot2)
library(tidyr)
library(dplR)
library(forecast)
library(imputeTS)
library(dlm)
library(itsmr)
library(MASS)

#*******************************************************************************
# Preparing the time series

raw_rwl <- readRDS("data/rwl_900+.Rds")

depth_function <- function(x){
    sum(!is.na(x))
}

rwl_mean <- apply(raw_rwl,1,mean,na.rm=T)
rwl_sd <- apply(raw_rwl,1,sd,na.rm=T)
rwl_depth <- apply(rwl_df, 1, depth_function)

rwl_ts <- ts(rwl_mean,end=2017)
rwl_mean_ts <- ts(rwl_mean,end=2017)
rwl_sd_ts <- ts(rwl_sd,end=2017)
rwl_depth_ts <- ts(rwl_depth,end=2017)

year <- seq(start(rwl_ts)[1],end(rwl_ts)[1])
seriesDepth <- as.data.frame(cbind(year,rwl_depth))
seriesMean <- as.data.frame(cbind(year,rwl_mean))

seriesDepthPlot <- ggplot(seriesDepth, aes(x=year,y=rwl_depth)) +
    geom_line() +
    ylab("Sample depth [#]") +
    xlab("Time") +
    theme(aspect.ratio = 0.618)

seriesMeanPlot <- ggplot(seriesMean, aes(x=year,y=rwl_mean)) +
    geom_line() +
    ylab("Mean ringwidth [1/100 mm]") +
    xlab("Time") +
    theme(aspect.ratio = 0.618)

# Showing off heteroskedasticity
heterosk <- lm(seriesMean$rwl_mean~seriesMean$year)
seriesHeterosk <- as.data.frame(cbind(year,heterosk$residuals))
seriesHeteroskedasticityPlot <- ggplot(seriesHeterosk, aes(x=year,y=V2)) +
    geom_line() +
    ylab("Residuals [1/100 mm]") +
    xlab("Time") +
    theme(aspect.ratio = 0.618)

spruce_window <- window(rwl_ts,start=1400, end=1800)
rwl_depth_window <- window(rwl_depth_ts,start=1400, end=1800)
rwl_mean_window <- window(rwl_mean_ts,start=1400, end=1800)
rwl_sd_window <- window(rwl_sd_ts,start=1400, end=1800)
t <- 1400:1800

data_to_plot <- as.data.frame(cbind(Time=t, ringwidth=spruce_window,model=rep("Observed",times=length(t))))

# Literature
# New Zealand Journal of Ecology (1990) 13: 9-15 (https://newzealandecology.org/nzje/1872)


# Approaches 1 + 2 are derived from literature
# Approaches 3 + 4 are more "intuitive".

# Approach 4 is the only one which displays exponential decrease of the ACF.


# We don't know which approach would be the way to go forward. -> Matthieu?


################################################################################
################################################################################
# Making series stationary: First approach
################################################################################
################################################################################


# Methods to estimate trend
################################################################################


# Method 1: Proposed in Woollons and Norton
#-------------------------------------------------------------------------------
# Warren (1980) (cited in Woollons and Norton, 1990) propose
# to estimate a trend of the form:
# Y = alpha*t^(beta)*exp(delta*t)
# This approach is nice because, if we
# take the ln on both sides:
# ln(Y) = ln(alpha) + beta*ln(t) + delta*t
# we get something linear in ln(t)

log_spruce_window <- log(spruce_window)
log_t <- log(t)
d1 <- cbind(log_spruce_window, log_t, t)

m1 <- lm(log_spruce_window~log_t+t, data=d1)
summary(m1)

# Fitted values
app_1_trend_1_fit<- exp(m1$fitted.values)

# ggplot
d_ggplot_1 <- cbind(Time=t, y1=spruce_window, y2=app_1_trend_1_fit)
d_ggplot_1 <- as.data.frame(d_ggplot_1)

pplot <- ggplot(d_ggplot_1, aes(x=Time))
pplot + geom_line(aes(y=y1), color="blue") +
    geom_line(aes(y=y2), color="black") +
    ggtitle("Estimated mean by the method of Warren (1980)")  +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab("Time") + ylab("Tree ring width (1/100 mm)")



# Method 2: Polynomial of order 2
##------------------------------------------------------------------------------
t2 <- t^2
d2 <- cbind(d1,t2)
d2 <- cbind(spruce_window, d2)
colnames(d2) <- c("y", "ln_y", "ln_t", "t", "t2")
head(d2)

m2 <- lm(y ~ t+t2,data=d2)
summary(m2)

app_1_trend_2_fit <- m2$fitted.values

# ggplot
d_ggplot_2 <- cbind(Time=t, y1=spruce_window, y2=app_1_trend_1_fit, y3=app_1_trend_2_fit)
d_ggplot_2 <- as.data.frame(d_ggplot_2)

data_p_order2 <- as.data.frame(cbind(Time=t, y=app_1_trend_2_fit,model=rep("Trend Polynomial order 2",times=length(t))))

d_ggplot_2_gathered <- d_ggplot_2 %>% gather(y, values, y2:y3)
colnames(d_ggplot_2_gathered) <- c("Time", "Width", "Est_Mean_Method", "Values")
str(d_ggplot_2_gathered)

pplot <- ggplot(d_ggplot_2_gathered, aes(x=Time))
pplot + geom_line(aes(y=Values, group=Est_Mean_Method, color=factor(Est_Mean_Method, labels=c("Warren (1980)", "Polynomial of order 2")),
                      linetype=factor(Est_Mean_Method, labels=c("Warren (1980)", "Polynomial of order 2")))) +
    geom_line(aes(y=Width), color="blue") +
    scale_color_manual(values=c("black", "green")) +
    ggtitle("Estimated means by the method of Warren (1980) \n and by a polynomial of order 2")  +
    theme(plot.title = element_text(hjust=0.5), legend.position = "top") +
    xlab("Time") + ylab("Tree ring width (1/100 mm)") +
    labs(color = "Methods") +
    labs(linetype= "Methods")


# Variance Stabilisation
################################################################################

# Residual time series
#-------------------------------------------------------------------------------

m1_res <- spruce_window - app_1_trend_1_fit
m2_res <- spruce_window - app_1_trend_2_fit

# ggplot
d_ggplot_3 <- cbind(Time=t, y1=m1_res, y2=m2_res)
d_ggplot_3 <- as.data.frame(d_ggplot_3)
str(d_ggplot_3)

d_ggplot_3_gathered <- d_ggplot_3 %>% gather(y, values, y1:y2)
str(d_ggplot_3_gathered)

# Legend manipulations
# https://rpubs.com/hughes/10012
pplot <- ggplot(d_ggplot_3_gathered, aes(x=Time))
pplot + geom_line(aes(y=values, group=y, color=factor(y, labels=c("Warren (1980)", "Polynomial of order 2")),
                      linetype=factor(y, labels=c("Warren (1980)", "Polynomial of order 2")))) +
    scale_color_manual(values=c("black", "green")) +
    ggtitle('Residual time series')  +
    theme(plot.title = element_text(hjust=0.5), legend.position = "top") +
    xlab("Time") + ylab("Tree ring width (1/100 mm)") +
    labs(color = "Methods") +
    labs(linetype="Methods")



# Residual transformations
#-------------------------------------------------------------------------------


# Method 1: Proposed transformation in Wollons and Norton (1990) for both methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1_i <- m1_res/app_1_trend_1_fit # Warren (1980)
m2_i <- m2$residuals/m2$fitted.values # Polynomial of order 2

temp <- cbind(Time=t, y=m1_i,model=rep("Woollons",times=length(t)))
data_woollons <- as.data.frame(cbind(Time=t, y=m1_i,model=rep("Woollons",times=length(t))))
data_p_order2 <- as.data.frame(cbind(Time=t, y=m2_i,model=rep("Polynomial order 2",times=length(t))))

# data for ggplot
d_ggplot_4 <- cbind(Time=t, y1=m1_i, y2=m2_i)
d_ggplot_4 <- as.data.frame(d_ggplot_4)
d_ggplot_4_gathered <- d_ggplot_4 %>% gather(y, values, y1:y2)



# Acf and Pacf
# For the acf and pacf plots the transformed residuals calculated with the method of Warren (1980)
# are used.
ggAcf(m1_i, main="ACF plot of the transformed residuals after having detrended the series with the method \n of Warren (1980)") # bold font does not work
ggPacf(m1_i, main="PACF plot of the transformed residuals after having detrended the series with the method \n of Warren (1980)")



# Method 2: Stabilisation by using SD of yearly observations (raw data)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Stabilisation: (sample mean - \mu)/sd(sample mean)
# \mu was estimated using Woollons and Norton

stab_sd_ts <- (rwl_mean_window - app_1_trend_1_fit)/(sqrt(rwl_depth_window)^(-1)*rwl_sd_window)
plot(stab_sd_ts)

d_ggplot_5 <- cbind(Time=t, y=stab_sd_ts)
d_ggplot_5 <- as.data.frame(d_ggplot_5)

data_woollons_sd_stab <- as.data.frame(cbind(Time=t, y=stab_sd_ts,model=rep("Woollons SD stab.",times=length(t))))
data_woollons_sd_stab$Time <- as.numeric(levels(data_woollons_sd_stab$Time)[data_woollons_sd_stab$Time])
data_woollons_sd_stab$y <- as.numeric(levels(data_woollons_sd_stab$y)[data_woollons_sd_stab$y])


################################################################################
################################################################################
# Second Approach: Power transformation of the residuals (by hand)
################################################################################
################################################################################

eps <- 1/10000 # To prevent log from getting -inf
R <- spruce_window
local_mean <- (R[1:(length(R)-1)] + R[2:length(R)])/2
local_sd <- abs(diff(R)) + eps

log_local_mean <- log(local_mean)
log_local_sd <- log(local_sd)

powt_df <- as.data.frame(cbind(log_S=log_local_sd, log_M=log_local_mean))
str(powt_df)

# linear model
m <- lm(log_S ~ log_M, data=powt_df)
summary(m)

b <- m$coef[2]

R_transformed <- R^(1-b)
data_powertrans <- as.data.frame(cbind(Time=t, y=R_transformed,method=rep("Power transformed",times=length(t))))


plot(y=R_transformed, x=1:length(R_transformed), type="l")

# Time series still has a trend
R_transformed_df <- data.frame(Width_trans=R_transformed, Time=t)
m <- lm(Width_trans~., data=R_transformed_df)
summary(m)

R_transformed_centered <- m$residuals

data_powertrans_centered <- as.data.frame(cbind(Time=t, y=R_transformed_centered,method=rep("Power trans. centered",times=length(t))))

acf(R_transformed_centered)
pacf(R_transformed_centered) # => AR(2)

arima(x=R_transformed_centered, c(2,0,0), order=c(0,0,0),method="ML")

# ggplot
# transformed time series
R_transformed
d_ggplot_6<- cbind(Time=t, y1=R_transformed, y2=R_transformed_centered)
d_ggplot_6<- as.data.frame(d_ggplot_6)
str(d_ggplot_5)

d_ggplot_5_gathered <- d_ggplot_6%>% gather(y, Values, y1:y2)

pplot <- ggplot(d_ggplot_5_gathered, aes(x=Time))
pplot + geom_line(aes(y=Values, group=y,
                      color=factor(y, labels=c("After power transformation", "After power transformation and linear trend")))) +
    scale_color_manual(values=c("blue", "green")) +
    ggtitle('Residual time series')  +
    theme(plot.title = element_text(hjust=0.5), legend.position = "top") +
    xlab("Time") + ylab("Tree ring width (1/100 mm)") +
    labs(color = "Methods") +
    labs(linetype="Methods")

ggAcf(R_transformed_centered, main="ACF plot after a power transformation and subtracting a linear trend")
ggPacf(R_transformed_centered, main="PACF plot after a power transformation and subtracting a linear trend")

################################################################################
################################################################################
# Third approach: Box-Cox transformation
################################################################################
################################################################################

# Box-Cox-Transform:
lambdas <- boxcox(spruce_window~t)
l <- lambdas$x[which.max(lambdas$y)] # this is the MLE lambda to transform data

spruce_window_boxcox <- (spruce_window^l-1)/l-mean((spruce_window^l-1)/l) # Box-Cox transformation
plotc(spruce_window_boxcox)

data_boxcox <- as.data.frame(cbind(Time=t, y=spruce_window_boxcox,method=rep("Box-Cox transform",times=length(t))))

acf(spruce_window_boxcox)
pacf(spruce_window_boxcox)


################################################################################
################################################################################
# Fourth approach: Log transform, polynomials
################################################################################
################################################################################

spruce_window_stab <- spruce_window/rwl_sd_window
spruce_window_log <- log(spruce_window_stab)

t2 <- t^2

order1_model <- lm(spruce_window_log ~ t)
summary(order1_model)

m2 <- lm(spruce_window_log ~ t+t2)
summary(m2)
# Nope, order 1 seems okay.
plotc(order1_model$residuals)

data_log_order1_varstab <- as.data.frame(cbind(Time=t, y=order1_model$residuals,model=rep("Log trans., order 1",times=length(t))))
data_log_order1_varstab$Time <- as.numeric(levels(data_log_order1_varstab$Time)[data_log_order1_varstab$Time])
data_log_order1_varstab$y <- as.numeric(levels(data_log_order1_varstab$y)[data_log_order1_varstab$y])



#*******************************************************************************
# Overview plots

# woollons
stationarity_woollons_plot <- ggplot(data_woollons_sd_stab,aes(x=Time,y=y)) +
    geom_line(size=0.5) +
    ylab("")

# power transformation & box-cox
data_powerbox <- rbind(data_powertrans_centered,
                       data_boxcox)
data_powerbox$Time <- as.numeric(levels(data_powerbox$Time)[data_powerbox$Time])
data_powerbox$y <- as.numeric(levels(data_powerbox$y)[data_powerbox$y])
stationarity_powerbox_plot <- ggplot(data_powerbox,aes(x=Time,y=y,colour=method)) +
    geom_line(size=0.5) +
    ylab("") +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

stationaritylog_order1_plot <- ggplot(data_log_order1_varstab,aes(x=Time,y=y)) +
    geom_line(size=0.5) +
    ylab("Scaled, log tranformed values [-]")
