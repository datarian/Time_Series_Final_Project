# getting window of the two series from 1332 to 1800:

library(ggplot2)
library(tidyr)
library(dplR)
library(forecast)
library(imputeTS)
library(dlm)

#spruce_window <- window(spruce_sup_900_ts,start=1400, end=1800)
raw_rwl <- readRDS("data/rwl_900+.Rds")
rwl_mean <- apply(raw_rwl,1,mean,na.rm=T)
rwl_sd <- apply(raw_rwl,1,sd,na.rm=T)
rwl_ts <- ts(rwl_mean,end=2017)
spruce_window <- window(rwl_ts,start=1400, end=1800)
t <- 1400:1800

# Literature
# New Zealand Journal of Ecology (1990) 13: 9-15 (https://newzealandecology.org/nzje/1872)


# Method proposed in Woollons and Norton
#########################################
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
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 124.270809  18.734549   6.633 1.07e-10 ***
#   log_t       -18.379202   2.940113  -6.251 1.05e-09 ***
#   t             0.010174   0.001844   5.516 6.25e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.1389 on 398 degrees of freedom
# Multiple R-squared:  0.5786,	Adjusted R-squared:  0.5765
# F-statistic: 273.2 on 2 and 398 DF,  p-value: < 2.2e-16

# Fitted values
m1_fitted<- exp(m1$fitted.values)
# plot(t, spruce_window, type="l")
# lines(t, m1_fitted, col="blue")

# ggplot
d_ggplot_1 <- cbind(Time=t, y1=spruce_window, y2=m1_fitted)
d_ggplot_1 <- as.data.frame(d_ggplot_1)
str(d_ggplot_1)
pplot <- ggplot(d_ggplot_1, aes(x=Time))
pplot + geom_line(aes(y=y1), color="blue") +
  geom_line(aes(y=y2), color="black") +
  ggtitle("Estimated mean by the method of Warren (1980)")  +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Time") + ylab("Tree ring width (mm)")

######################################################

# Polynomial of order 2
##########################
t2 <- t^2
d2 <- cbind(d1,t2)
d2 <- cbind(spruce_window, d2)
colnames(d2) <- c("y", "ln_y", "ln_t", "t", "t2")
head(d2)
# y     ln_y     ln_t    t      t2
# [1,] 132 4.882802 7.244228 1400 1960000
# [2,] 135 4.905275 7.244942 1401 1962801
# [3,] 183 5.209486 7.245655 1402 1965604
# [4,] 180 5.192957 7.246368 1403 1968409
# [5,] 157 5.056246 7.247081 1404 1971216
# [6,] 164 5.099866 7.247793 1405 1974025

m2 <- lm(y ~ t+t2,data=d2)
summary(m2)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  2.275e+03  2.348e+02   9.689  < 2e-16 ***
#   t           -2.449e+00  2.945e-01  -8.314 1.49e-15 ***
#   t2           6.975e-04  9.199e-05   7.582 2.41e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 22.08 on 398 degrees of freedom
# Multiple R-squared:  0.5906,	Adjusted R-squared:  0.5886
# F-statistic: 287.1 on 2 and 398 DF,  p-value: < 2.2e-16

m2_fitted <- m2$fitted.values
# plot(t, spruce_window, type="l")
# lines(t, m2_fitted, col="red")
# lines(t, m1_fitted, col="blue")

# ggplot
d_ggplot_2 <- cbind(Time=t, y1=spruce_window, y2=m1_fitted, y3=m2_fitted)
d_ggplot_2 <- as.data.frame(d_ggplot_2)
d_ggplot_2_gathered <- d_ggplot_2 %>% gather(y, values, y2:y3)
tail(d_ggplot_2_gathered)
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


# Without groups
# pplot <- ggplot(d_ggplot_2, aes(x=Time))
# pplot + geom_line(aes(y=y1), color="blue") +
#   geom_line(aes(y=y2), color="black") +
#   geom_line(aes(y=y3), color="green") +
#   ggtitle('Estimated means by the method of Warren (1980) in black \n and by a polynomial of order 2 in green')  +
#   theme(plot.title = element_text(hjust=0.5)) +
#   xlab("Time") + ylab("Tree ring width (mm)")

# Residual time series
########################

m1_res <- spruce_window - m1_fitted
m2_res <- spruce_window - m2_fitted

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
  #scale_linetype_discrete(name = "Fancy Title")
  #guide_legend(title="Method")


# Residual transformations
##########################

# Proposed transformation in Wollons and Norton (1990) for both methods
# ---------------------------------------------------------------------
m1_i <- m1_res/m1_fitted # Warren (1980)
m2_i <- m2$residuals/m2$fitted.values # Polynomial of order 2

# ymin <- min(min(m1_i), min(m2_i))
# ymax <- max(max(m1_i), max(m2_i))
# plot(t, m1_i, col="blue", type="l", ylim=c(ymin, ymax))
# lines(t, m2_i, col="red", lty=1)

# ggplot
d_ggplot_4 <- cbind(Time=t, y1=m1_i, y2=m2_i)
d_ggplot_4 <- as.data.frame(d_ggplot_4)
d_ggplot_4_gathered <- d_ggplot_4 %>% gather(y, values, y1:y2)

pplot <- ggplot(d_ggplot_4_gathered, aes(x=Time))
pplot + geom_line(aes(y=values, group=y, color=factor(y, labels=c("Warren (1980)", "Polynomial of order 2")),
                      linetype=factor(y, labels=c("Warren (1980)", "Polynomial of order 2")))) +
  scale_color_manual(values=c("black", "green")) +
  ggtitle('Transformed residuals')  +
  theme(plot.title = element_text(hjust=0.5), legend.position = "top") +
  xlab("Time") + ylab("Tree ring width (1/100 mm)") +
  labs(color = "Methods") +
  labs(linetype= "Methods")


# Acf and Pacf
# For the acf and pacf plots the transformed residuals calculated with the method of Warren (1980)
# are used.
ggAcf(m1_i, main="ACF plot of the transformed residuals after having detrended the series with the method \n of Warren (1980)") # bold font does not work
ggPacf(m1_i, main="PACF plot of the transformed residuals after having detrended the series with the method \n of Warren (1980)")


# Power transformation of the residuals (by hand)
# -----------------------------------------------
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

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -8.2373     2.5306  -3.255  0.00123 **
#   log_M         2.0442     0.5055   4.044 6.31e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b <- m$coef[2]

R_transformed <- R^(1-b)

plot(y=R_transformed, x=1:length(R_transformed), type="l")

# Time series still has a trend
R_transformed_df <- data.frame(Width_trans=R_transformed, Time=t)
m <- lm(Width_trans~., data=R_transformed_df)
summary(m)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) -6.083e-03  5.579e-04  -10.90   <2e-16 ***
#   Time         7.260e-06  3.478e-07   20.88   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.0008062 on 399 degrees of freedom
# Multiple R-squared:  0.522,	Adjusted R-squared:  0.5208
# F-statistic: 435.7 on 1 and 399 DF,  p-value: < 2.2e-16

R_transformed_demeaned <- m$residuals
# plot(y=R_transformed_demeaned, x=1:length(R_transformed_demeaned), type="l")
#
# acf(R_transformed_demeaned)
# pacf(R_transformed_demeaned) # => AR(2)

arima(x=R_transformed_demeaned, c(1,0,0), method="ML")
# Coefficients:
#   ar1  intercept
# 0.6255      0e+00
# s.e.  0.0393      1e-04
#
# sigma^2 estimated as 3.967e-07:  log likelihood = 2386.12,  aic = -4766.25

# ggplot
# transformed time series
R_transformed
d_ggplot_5 <- cbind(Time=t, y1=R_transformed, y2=R_transformed_demeaned)
d_ggplot_5 <- as.data.frame(d_ggplot_5)
str(d_ggplot_5)

d_ggplot_5_gathered <- d_ggplot_5 %>% gather(y, Values, y1:y2)

pplot <- ggplot(d_ggplot_5_gathered, aes(x=Time))
pplot + geom_line(aes(y=Values, group=y, color=factor(y, labels=c("After power transformation", "After power transformation and linear trend")))) +
  scale_color_manual(values=c("blue", "green")) +
  ggtitle('Residual time series')  +
  theme(plot.title = element_text(hjust=0.5), legend.position = "top") +
  xlab("Time") + ylab("Tree ring width (1/100 mm)") +
  labs(color = "Methods") +
  labs(linetype="Methods")

ggAcf(R_transformed_demeaned, main="ACF plot after a power transformation and subtracting a linear trend")
ggPacf(R_transformed_demeaned, main="PACF plot after a power transformation and subtracting a linear trend")

# Missing values
################
R_transformed_demeaned_nas <- R_transformed_demeaned

# Create missing values
R_transformed_demeaned_nas[120:150] <- NA

# imputTS package
smoothed_values_imputeTS <- na.kalman(x=R_transformed_demeaned_nas, model="auto.arima")

plot(y=R_transformed_demeaned, x=1:length(R_transformed_demeaned), type="l")
lines(y=y, x=1:length(R_transformed_demeaned), col="blue")

# dlm model
acf(R_transformed_demeaned_nas, na.action=na.pass)
pacf(R_transformed_demeaned_nas, na.action=na.pass) # => AR(2), with coeff

# AR(2) in state space form (https://robjhyndman.com/talks/ABS3.pdf):
# Let x_t = [y_t, y_t-1]^T, and w_t = [e_t,0]
# y_t = [1 0]x_t + v_t
# x_t = |phi_1 phi_2| x_{t-1} + w_t
#       |1        0 |

# In notation of the dlm package:
# FF = [1 0]
# GG = |phi_1 phi_2|
#      |1      0   |
# An R Package for Dynamic Linear Models, p.9
# parameters x:
# x[1] = phi_1
# x[2] = phi_2
# x[3] = sqrt(var(v_t))
# x[4] = sqrt(var(w_t))

build <- function(x){
  #x <- c(0.6, 0.25, 0.01, 0.01, 0.03)
  FF <- matrix(c(1,0), nrow=1)
  GG <- matrix(c(x[1], x[2], 1, 0), byrow=T, nrow=2)
  V = matrix(x[3]^2,nrow=1)
  W = crossprod(diag(c(x[4], x[5])))

  mod <- dlm(FF=FF, GG=GG, V=V, W=W, m0=matrix(c(0,0), ncol=1), C0=diag(c(10^7, 10^7)))
  return(mod)
}

# This gives almost the same as na.kalman
initial_params <- c(0.6, 0.25, 0.01, 0.01, 0.05) # maybe we have to come up with better initial values for the standard deviations
build(initial_params)

y <- as.numeric(R_transformed_demeaned_nas)
str(y)

fit2 <- dlmMLE(y, parm=initial_params, build=build, control = list(maxit = 10000))
fit2$par # 0.60 0.25 0.01 0.01 0.01
# [1]  5.998274e-01  2.402065e-01 -3.767355e-05  2.182666e-04 -2.531772e-03

str(fit2)

#filtered_values <- dlmFilter(R_transformed_demeaned_nas, build(fit2$par))$f
smoothed_states <- dlmSmooth(R_transformed_demeaned_nas, build(fit2$par))$s

smoothed_values <- t(matrix(c(1,0), nrow=1) %*% t(smoothed_states))

# Visual verification
plot(y=R_transformed_demeaned, x=1:length(R_transformed_demeaned), type="l")
lines(y=smoothed_values[-1], x=1:length(R_transformed_demeaned), col="red")
lines(y=smoothed_values_imputeTS, x=1:length(R_transformed_demeaned), col="blue", lty=2)

