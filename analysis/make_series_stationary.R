# getting window of the two series from 1332 to 1800:

library(ggplot2)
library(tidyr)
library(dplR)
library(forecast)


spruce_window <- window(spruce_sup_900_ts,start=1400, end=1800)
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
  xlab("Time") + ylab("Tree ring width (mm)") +
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

# Residual time series
# --------------------
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
  xlab("Time") + ylab("Tree ring width (mm)") +
  labs(color = "Methods") +
  labs(linetype="Methods")
  #scale_linetype_discrete(name = "Fancy Title")
  #guide_legend(title="Method")


# Residual transformations
##########################

# Proposed transformation in Wollons and Norton (1990) for both methods
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
  xlab("Time") + ylab("Tree ring width (mm)") +
  labs(color = "Methods") +
  labs(linetype= "Methods")


# Acf and Pacf
################

# For the acf and pacf plots the transformed residuals calculated with the method of Warren (1980) 
# are used.
ggAcf(m1_i, main="ACF plot of the transformed residuals after having detrended the series with the method \n of Warren (1980)") # bold font does not work
ggPacf(m1_i, main="PACF plot of the transformed residuals after having detrended the series with the method \n of Warren (1980)")



