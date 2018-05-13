################################################################################

# Rk: the residuals are taken from the fourth approach

################################################################################

#library(MASS)
# theta <- -0.8
# res <-mvrnorm(1, mu=c(0,0), Sigma=matrix(c(1,theta,theta,theta^2),2,2))
# b <- res[1]
# c <- res[2]
# c/b
#
# K <- matrix(c(1,theta,0,0), 2,2)
# W <- K%*%t(K)
# res <-mvrnorm(1, mu=c(0,0), Sigma=W)
# b <- res[1]
# c <- res[2]
# c/b


rm(list = ls())

library(forecast)
library(matrixcalc)

source("analysis/prepare_individual_series.R")
source("analysis//make_series_stationary.R")

# Literature
# ------------------------------------------------------------------------------
# https://sites.ualberta.ca/~sfossati/e509/files/other/dlm_ex.R
# https://robjhyndman.com/talks/ABS3.pdf

# Original transformed time series
#-------------------------------------------------------------------------------

ts_res <- m$residuals
mean(ts_res) 4.697421e-18

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Imuputations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# We formulate an AR(1), AR(2) and an ARMA(1,1) as state space models. Imputations
# are done with the Kalman smoother. We also check the results with the Kalman
# smoother of th imputeTS package

# set.seed(123)
# # We divide the time series in 40 blocks of 9 data points. The blocks are separated by one
# # value. This ensures that never two missing values occur in sequence
# # | 1, 2,..., 9| 10 | 11, 13, ..., 29 | 30 | 31, ..., 39| ... | 390 | 391, ..., 399| 400|
# #impute_idx <- sample(1:401, 25)
# impute_idx <- sample(1:9, 40, replace=T) + seq(0, 390, by=10)
#
# # Transform time series
# # ------------------------------------------------------------------------------
# t <- 1400:1800
# t <- t[-impute_idx]
# spruce_window_stab <- spruce_window[-impute_idx]/rwl_sd_window[-impute_idx]
# spruce_window_log <- log(spruce_window_stab)
#
# m <- lm(spruce_window_log ~ t)
# summary(m)
#
# ts_res <- m$residuals
# mean(ts_res) # 5.276676e-18

# ------------------------------------------------------------------------------

# We also need define the process of rwl_sd_window. For simplicity we use the
# original times series (no NAs)
plot(y=rwl_sd_window, x=1:length(rwl_sd_window), type="l")
acf(rwl_sd_window)
pacf(rwl_sd_window) # AR(2)

mean_sd <- mean(rwl_sd_window)
rwl_sd_window_centered <- as.numeric(rwl_sd_window - mean_sd)

# ------------------------------------------------------------------------------

# Best models according to the auto.fit function
################################################################################
auto.arima(ts_res, d=0) # => ARIMA(1,0,1)
auto.arima(rwl_sd_window_centered, d=0) # => ARIMA(1,0,1)

# Models
################################################################################

# AR(1) state space model
#-------------------------------------------------------------------------------
# AR(1):
# x_t = phi*x_{t-1} + eta_t, eta_t ~ N(0, sigma^2_{eta})

# y_t = x_t
# x_t = x_{t+1} + w_t, w_t ~ N(0,W_t)
# x_0 ~ N(m_0, C_0), with m_0 = 0


# Rk: sigma^2_x = sigma^2_{eta}/(1 -phi^2)

buildAR1 <- function(x){

    phi <- x[1] # parameter of the AR(1) process
    sd <- x[2]

    FF <- 1
    V <- 0
    GG <- phi
    W <- sd^2

    #C0 <- sd^2/(1-phi^2)

    mAR1 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=0, C0=1e4)
    return(mAR1)
}


# AR(2) as state space model
#-------------------------------------------------------------------------------
# AR(2) in state space form (https://robjhyndman.com/talks/ABS3.pdf):
# Let x_t = [y_t, y_t-1]^T, and w_t = [e_t,0] and v_t = 0
# y_t = [1 0]x_t + v_t
# x_t = |phi_1 phi_2| x_{t-1} + w_t
#       |1        0 |

# In notation of the dlm package:
# FF = [1 0]
# GG = |phi_1 phi_2|
#      |1      0   |

buildAR2 <- function(x){

    phi1 <- x[1]
    phi2 <- x[2]
    sd <- x[3]

    FF <- matrix(c(1,0), nrow=1)
    GG <- matrix(c(phi1, phi2, 1, 0), byrow=T, nrow=2)
    V <- 0 #1e-5 # Due to convergence issues we define V to be slightly different from 0
    W <- diag(c(sd^2, 0))

    mAR2 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=matrix(c(0,0), ncol=1), C0=diag(c(1e4, 1e4)))
    return(mAR2)
}

# ARMA(1,1) as state space model
#-------------------------------------------------------------------------------
# y_t = phi*y_{t-1} + theta*e_{t-1}

# Let:
# x_t = |y_t        |
#       |theta * e_t|

# and:
# w_t = |e_t      |
#       |theta*e_t|

# Model:
# y_t = [1 0]*x_t
# x_t = |phi 1|*x_{t-1} + w_t
#       |0   0|


buildARMA11 <- function(x){
    phi <- x[1]
    theta <- x[2]
    # if (abs(theta) > 1) {
    #     #theta <- 1/theta
    #     print("yes")
    # }

    sd <- x[3]
    FF <- matrix(c(1,0), nrow=1)
    GG <- matrix(c(phi,0,1,0), 2,2)

    V <- 0

    #theta <- -1.3531236
    #sd <- 0.09245244
    K_W <- matrix(c(1,theta,0,0), 2,2) # K_W SHOULD BE LIKE THIS, BUT DLM HATES IT, BECAUSE IT APPLIES A CHOLESKY DECOMPOSITION
    K_W <- K_W + matrix(c(1e-5,1e-5,1e-5,1e-5),2,2)
    W <- (K_W%*%t(K_W))%*%diag(c(sd^2, sd^2))
    # print(is.positive.semi.definite(W))
    #print(chol(W))

    mARMA11 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=matrix(c(0,0), ncol=1), C0=diag(c(1e4, 1e4)))
    return(mARMA11)
}


# We also need to estimate the variance. We model it as an AR(2) process

# Ex. start values
################################################################################
# Rk: ts_res has already mean 0

# AR(1)
#-------------------------------------------------------------------------------
AR1 <- arima(ts_res, order=c(1,0,0), include.mean = F)
phi_AR1_start <- as.numeric(AR1$coef) # parameter 1 to be optimized
#residuals_sd <- sqrt(400/401*var(ar1$residuals)) # parameter2 to be optimized
sigma_sd <- sqrt(as.numeric(AR1$sigma2))

################################################################################
# Simulation study 1: Standard deviation is estimated as an AR(2) process
################################################################################

# We compare the AR(1), AR(2) and ARMA(1,1) with respect to how good they are
# at imputing missing values.

# One value to impute
#-------------------------------------------------------------------------------
set.seed(123)
# We divide the time series in 40 blocks of 9 data points. The blocks are separated by one
# value. This ensures that never two missing values occur in sequence
# | 1, 2,..., 9| 10 | 11, 13, ..., 29 | 30 | 31, ..., 39| ... | 390 | 391, ..., 399| 400|
#impute_idx <- sample(1:401, 25)
impute_idx <- sample(1:9, 40, replace=T) + seq(0, 390, by=10)

imputed_widths_untransformed_1val <- matrix(NA, nrow=3, ncol=length(impute_idx))
imputed_widths_transformed_1val <- imputed_widths_untransformed_1val
imputed_sd_1val <- rep(NA, length(impute_idx))

convergence_mat <- matrix(NA, nrow=4, ncol=length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- ts_res
    y_sd <- rwl_sd_window_centered
    y[i] <- NA
    y_sd[i] <- NA

    # Imputation
    ############################################################################

    # Tree ring width
    #---------------------------------------------------------------------------
    # AR(1)
    #...........................................................................
    AR1 <-arima(y, order=c(1,0,0), include.mean = F) # We assume that the series has still mean 0
    phi <- as.numeric(AR1$coef)
    sigma <- sqrt(as.numeric(AR1$sigma2))

    initial_params <- c(phi, sigma)
    fit <- try(dlmMLE(y, parm=initial_params, build=buildAR1, control = list(maxit = 1000), method="BFGS"))
    if(class(fit)=="try-error"){
        convergence_mat[2,j] <- 1
    } else {
        smoothed_state <- dlmSmooth(y, buildAR1(fit$par))$s[-1]
        imputed_widths_transformed_1val[1,j] <- smoothed_state[i] # F = 1
        convergence_mat[1,j] <- fit$convergence
    }

    # AR(2)
    #...........................................................................
    AR2 <-arima(y, order=c(2,0,0), include.mean = F) # We assume that the series has still mean 0
    phi1 <- as.numeric(AR2$coef[1])
    phi2 <- as.numeric(AR2$coef[2])
    sigma <- sqrt(as.numeric(AR2$sigma2))

    initial_params <- c(phi1, phi2, sigma)
    fit <- try(dlmMLE(y, parm=initial_params, build=buildAR2, control = list(maxit = 1000), method="BFGS"))
    if(class(fit)=="try-error"){
        convergence_mat[2,j] <- 1
    } else {
        smoothed_states <- dlmSmooth(y, buildAR2(fit$par))$s
        smoothed_state1 <- smoothed_states[-1,1][i]
        smoothed_state2 <- smoothed_states[-1,2][i]
        imputed_widths_transformed_1val[2,j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat[2,j] <- fit$convergence
    }

    # ARMA(1,1)
    #...........................................................................
    ARMA11 <-arima(y, order=c(1,0,1), include.mean = F) # We assume that the series has still mean 0
    phi <- as.numeric(ARMA11$coef[1])
    theta <- as.numeric(ARMA11$coef[2])
    sigma <- sqrt(as.numeric(ARMA11$sigma2))

    initial_params <- c(phi, theta, sigma)
    fit <- try(dlmMLE(y, parm=initial_params, build=buildARMA11, control = list(maxit = 1000), method="BFGS"))
    if(class(fit)=="try-error"){
        convergence_mat[3,j] <- 1
    } else {
        smoothed_states <- dlmSmooth(y, buildARMA11(fit$par))$s
        smoothed_state1 <- smoothed_states[-1,1][i]
        smoothed_state2 <- smoothed_states[-1,2][i]
        imputed_widths_transformed_1val[3,j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat[3,j] <- fit$convergence
    }

    # Tree ring sd
    #---------------------------------------------------------------------------
    AR2 <-arima(y_sd, order=c(2,0,0), include.mean = F) # We assume that the series has still mean 0
    phi1 <- as.numeric(AR2$coef[1])
    phi2 <- as.numeric(AR2$coef[2])
    sigma <- sqrt(as.numeric(AR2$sigma2))

    initial_params <- c(phi1, phi2, sigma)
    fit <- try(dlmMLE(y_sd, parm=initial_params, build=buildAR2, control = list(maxit = 1000), method="BFGS"))
    if(class(fit)=="try-error"){
        convergence_mat[4,j] <- 1
    } else {
        smoothed_states <- dlmSmooth(y_sd, buildAR2(fit$par))$s
        smoothed_state1 <- smoothed_states[-1,1][i]
        smoothed_state2 <- smoothed_states[-1,2][i]
        imputed_sd_1val[j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat[4,j] <- fit$convergence
    }
    j <- j + 1
}


imputed_widths_transformed_1val
imputed_sd_1val <- imputed_sd_1val + mean_sd
convergence_mat # 0 means that optimization algorithm has converged


# Back transformation
# ------------------------------------------------------------------------------
t <- 1400:1800
#t2 <- t^2
coefs <- matrix(as.numeric(m$coef), nrow=2)
x <- cbind(rep(1, times=length(ts_res)), t)

log_val_mean <-x[impute_idx,]%*%coefs # alpha + beta*t

# ------------------------------------------------------------------------------

tmp <- imputed_widths_transformed_1val + matrix(log_val_mean, ncol=length(impute_idx), nrow=3,byrow=T)

# transformed time series: imputed vs. true
plot(y=imputed_widths_transformed_1val[1,], x= t[impute_idx], type="l", col="blue") # estimated
lines(y=ts_res[impute_idx], x=t[impute_idx]) # true

# add linear trend back
plot(y=exp(tmp[1,]), x= t[impute_idx], type="l", col="blue") # estimated
lines(y=exp(ts_res[impute_idx] + log_val_mean), x=t[impute_idx]) # true

# average sd used to calculate estimated tree widths (bad as expected)
plot(y=exp(tmp[1,])*mean_sd, x= t[impute_idx], type="l", col="blue")  # estimated
lines(y=exp(ts_res[impute_idx] + log_val_mean)*(rwl_sd_window_centered[impute_idx]+mean_sd), x=t[impute_idx]) # true

# estimated sd
plot(y=exp(tmp[1,])*(mean_sd+imputed_sd_1val-mean_sd), x= t[impute_idx], type="l", col="blue") # estimated, with estimated sd
lines(y=exp(ts_res[impute_idx] + log_val_mean)*(rwl_sd_window_centered[impute_idx]+mean_sd), x=t[impute_idx]) # true

# true sd
plot(y=exp(tmp[1,])*(rwl_sd_window_centered[impute_idx]+mean_sd), x= t[impute_idx], type="l", col="blue") # estimated, with true sd
lines(y=exp(ts_res[impute_idx] + log_val_mean)*(rwl_sd_window_centered[impute_idx]+mean_sd), x=t[impute_idx]) # true

# ------------------------------------------------------------------------------

imputed_widths_untransformed_1val <- exp(imputed_widths_transformed_1val +
                                         matrix(log_val_mean, ncol=length(impute_idx),
                                                nrow=3,byrow=T))*(imputed_sd_1val)

# Transformed width, i.e before back transformation
# -------------------------------------------------
imputed_widths_transformed_1val_all <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_transformed_1val_all[,impute_idx] <- imputed_widths_transformed_1val
imputed_widths_transformed_1val_all[,impute_idx-1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_transformed_1val_all[,impute_idx+1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx+1]

# Width after back transformation
# -------------------------------
imputed_widths_1val <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_1val[,impute_idx] <- imputed_widths_untransformed_1val
imputed_widths_1val[,impute_idx-1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_1val[,impute_idx+1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx+1]

# True values
# ----------
true_values_transformed_1val <- matrix(NA, nrow=1, ncol= length(ts_resl))
true_values_transformed_1val[-impute_idx] <- ts_res[-impute_idx]
true_values_untransformed_1val <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val[-impute_idx] <- spruce_window[-impute_idx]

# This contains the true values at the points where the imputatiions are done
true_values_transformed_1val_2 <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val_2[impute_idx] <- ts_res[impute_idx]
true_values_transformed_1val_2[impute_idx-1] <- ts_res[impute_idx-1]
true_values_transformed_1val_2[impute_idx+1] <- ts_res[impute_idx+1]

true_values_untransformed_1val_2 <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val_2[impute_idx] <- spruce_window[impute_idx]
true_values_untransformed_1val_2[impute_idx-1] <- spruce_window[impute_idx-1]
true_values_untransformed_1val_2[impute_idx+1] <- spruce_window[impute_idx+1]

# Plots: transformed time series
plot(y=true_values_transformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_transformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_transformed_1val_all[1,], x=1:length(ts_res), col="green")
lines(imputed_widths_transformed_1val_all[2,], x=1:length(ts_res), col="blue")
lines(imputed_widths_transformed_1val_all[3,], x=1:length(ts_res), col="red")

# Plots: Untransformed (i.e. initial) time series
plot(y=true_values_untransformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_untransformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_1val[1,], x=1:length(ts_res), col="green", lty=2)
lines(imputed_widths_1val[2,], x=1:length(ts_res), col="blue", lty=2)
lines(imputed_widths_1val[3,], x=1:length(ts_res), col="red", lty=2)

# MSE transformed
################################################################################
mse_transformed_1val <- rowSums((imputed_widths_transformed_1val -
                                     matrix(ts_res[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                     length(impute_idx)
# 0.005795514 0.005581524 0.005365752

# MSE back transformed (untransformed)
################################################################################
mse_untransformed_1val <- rowSums((imputed_widths_untransformed_1val -
                                    matrix(spruce_window[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                    length(impute_idx)
# 3764.375 3896.546 4129.142


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


################################################################################
# Simulation study 2: Standard deviation is estimated as an ARMA(1,1) process
################################################################################

# We compare the AR(1), AR(2) and ARMA(1,1) with respect to how good they are
# at imputing missing values.

# One value to impute
#-------------------------------------------------------------------------------
set.seed(123)
# We divide the time series in 40 blocks of 9 data points. The blocks are separated by one
# value. This ensures that never two missing values occur in sequence
# | 1, 2,..., 9| 10 | 11, 13, ..., 29 | 30 | 31, ..., 39| ... | 390 | 391, ..., 399| 400|
impute_idx <- sample(1:9, 40, replace=T) + seq(0, 390, by=10)
imputed_sd_1val <- rep(NA, length(impute_idx))

convergence <- rep(NA, times=length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- ts_res
    y_sd <- rwl_sd_window_centered
    y[i] <- NA
    y_sd[i] <- NA

    # Imputation
    ############################################################################

    # Tree ring sd
    #---------------------------------------------------------------------------
    ARMA11 <-arima(y_sd, order=c(1,0,1), include.mean = F) # We assume that the series has still mean 0
    phi <- as.numeric(ARMA11$coef[1])
    theta <- as.numeric(ARMA11$coef[2])
    sigma <- sqrt(as.numeric(ARMA11$sigma2))

    initial_params <- c(phi, theta, sigma)
    fit <- try(dlmMLE(y_sd, parm=initial_params, build=buildARMA11, control = list(maxit = 1000), method="BFGS"))
    if(class(fit)=="try-error"){
        convergence_mat[4,j] <- 1

    } else {
        smoothed_states <- dlmSmooth(y_sd, buildARMA11(fit$par))$s
        smoothed_state1 <- smoothed_states[-1,1][i]

        smoothed_state2 <- smoothed_states[-1,2][i]
        imputed_sd_1val[j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence[j] <-fit$convergence
    }

    j <- j + 1
}


imputed_widths_transformed_1val
imputed_sd_1val <- imputed_sd_1val + mean_sd
convergence # 0 means that optimization algorithm has converged


# Back transformation
# ------------------------------------------------------------------------------
t <- 1400:1800
#t2 <- t^2
coefs <- matrix(as.numeric(m$coef), nrow=2)
x <- cbind(rep(1, times=length(ts_res)), t)

log_val_mean <-x[impute_idx,]%*%coefs # alpha + beta*t

imputed_widths_untransformed_1val <- exp(imputed_widths_transformed_1val +
                                             matrix(log_val_mean, ncol=length(impute_idx),
                                                    nrow=3,byrow=T))*(imputed_sd_1val)

# Transformed width, i.e before back transformation
# -------------------------------------------------
imputed_widths_transformed_1val_all <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_transformed_1val_all[,impute_idx] <- imputed_widths_transformed_1val
imputed_widths_transformed_1val_all[,impute_idx-1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_transformed_1val_all[,impute_idx+1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx+1]

# Width after back transformation
# -------------------------------
imputed_widths_1val <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_1val[,impute_idx] <- imputed_widths_untransformed_1val
imputed_widths_1val[,impute_idx-1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_1val[,impute_idx+1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx+1]

# True values
# ----------
true_values_transformed_1val <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val[-impute_idx] <- ts_res[-impute_idx]
true_values_untransformed_1val <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val[-impute_idx] <- spruce_window[-impute_idx]

# This contains the true values at the points where the imputatiions are done
true_values_transformed_1val_2 <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val_2[impute_idx] <- ts_res[impute_idx]
true_values_transformed_1val_2[impute_idx-1] <- ts_res[impute_idx-1]
true_values_transformed_1val_2[impute_idx+1] <- ts_res[impute_idx+1]

true_values_untransformed_1val_2 <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val_2[impute_idx] <- spruce_window[impute_idx]
true_values_untransformed_1val_2[impute_idx-1] <- spruce_window[impute_idx-1]
true_values_untransformed_1val_2[impute_idx+1] <- spruce_window[impute_idx+1]

# Plots: transformed time series
plot(y=true_values_transformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_transformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_transformed_1val_all[1,], x=1:length(ts_res), col="green")
lines(imputed_widths_transformed_1val_all[2,], x=1:length(ts_res), col="blue")
lines(imputed_widths_transformed_1val_all[3,], x=1:length(ts_res), col="red")

# Plots: Untransformed (i.e. initial) time series
plot(y=true_values_untransformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_untransformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_1val[1,], x=1:length(ts_res), col="green", lty=2)
lines(imputed_widths_1val[2,], x=1:length(ts_res), col="blue", lty=2)
lines(imputed_widths_1val[3,], x=1:length(ts_res), col="red", lty=2)

# MSE transformed
################################################################################
mse_transformed_1val <- rowSums((imputed_widths_transformed_1val -
                                     matrix(ts_res[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                     length(impute_idx)
# 0.005795514 0.005581524 0.005365759

# MSE back transformed (untransformed)
################################################################################
mse_untransformed_1val <- rowSums((imputed_widths_untransformed_1val -
                                       matrix(spruce_window[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                       length(impute_idx)
# 3482.894 3670.342 3961.353



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


################################################################################
# Simulation study 3: Standard deviation is estimated as the mean of the values
# observed before and after the missing value
################################################################################

# We compare the AR(1), AR(2) and ARMA(1,1) with respect to how good they are
# at imputing missing values.

# One value to impute
#-------------------------------------------------------------------------------
set.seed(123)
# We divide the time series in 40 blocks of 9 data points. The blocks are separated by one
# value. This ensures that never two missing values occur in sequence
# | 1, 2,..., 9| 10 | 11, 13, ..., 29 | 30 | 31, ..., 39| ... | 390 | 391, ..., 399| 400|
impute_idx <- sample(1:9, 40, replace=T) + seq(0, 390, by=10)
imputed_sd_1val <- rep(NA, length(impute_idx))

convergence <- rep(NA, times=length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- ts_res
    y_sd <- rwl_sd_window_centered
    y[i] <- NA
    y_sd[i] <- NA

    # Imputation
    ############################################################################

    # Tree ring sd
    #---------------------------------------------------------------------------
    imputed_sd_1val[j] <- (y_sd[i-1] + y_sd[i+1])/2

    j <- j + 1
}


imputed_widths_transformed_1val
imputed_sd_1val <- imputed_sd_1val + mean_sd
convergence_mat # 0 means that optimization algorithm has converged


# Back transformation
# ------------------------------------------------------------------------------
t <- 1400:1800
#t2 <- t^2
coefs <- matrix(as.numeric(m$coef), nrow=2)
x <- cbind(rep(1, times=length(ts_res)), t)

log_val_mean <-x[impute_idx,]%*%coefs # alpha + beta*t

imputed_widths_untransformed_1val <- exp(imputed_widths_transformed_1val +
                                             matrix(log_val_mean, ncol=length(impute_idx),
                                                    nrow=3,byrow=T))*(imputed_sd_1val)

# Transformed width, i.e before back transformation
# -------------------------------------------------
imputed_widths_transformed_1val_all <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_transformed_1val_all[,impute_idx] <- imputed_widths_transformed_1val
imputed_widths_transformed_1val_all[,impute_idx-1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_transformed_1val_all[,impute_idx+1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx+1]

# Width after back transformation
# -------------------------------
imputed_widths_1val <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_1val[,impute_idx] <- imputed_widths_untransformed_1val
imputed_widths_1val[,impute_idx-1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_1val[,impute_idx+1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx+1]

# True values
# ----------
true_values_transformed_1val <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val[-impute_idx] <- ts_res[-impute_idx]
true_values_untransformed_1val <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val[-impute_idx] <- spruce_window[-impute_idx]

# This contains the true values at the points where the imputatiions are done
true_values_transformed_1val_2 <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val_2[impute_idx] <- ts_res[impute_idx]
true_values_transformed_1val_2[impute_idx-1] <- ts_res[impute_idx-1]
true_values_transformed_1val_2[impute_idx+1] <- ts_res[impute_idx+1]

true_values_untransformed_1val_2 <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val_2[impute_idx] <- spruce_window[impute_idx]
true_values_untransformed_1val_2[impute_idx-1] <- spruce_window[impute_idx-1]
true_values_untransformed_1val_2[impute_idx+1] <- spruce_window[impute_idx+1]

# Plots: transformed time series
plot(y=true_values_transformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_transformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_transformed_1val_all[1,], x=1:length(ts_res), col="green")
lines(imputed_widths_transformed_1val_all[2,], x=1:length(ts_res), col="blue")
lines(imputed_widths_transformed_1val_all[3,], x=1:length(ts_res), col="red")

# Plots: Untransformed (i.e. initial) time series
plot(y=true_values_untransformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_untransformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_1val[1,], x=1:length(ts_res), col="green", lty=2)
lines(imputed_widths_1val[2,], x=1:length(ts_res), col="blue", lty=2)
lines(imputed_widths_1val[3,], x=1:length(ts_res), col="red", lty=2)

# MSE transformed
################################################################################
mse_transformed_1val <- rowSums((imputed_widths_transformed_1val -
                                     matrix(ts_res[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                     length(impute_idx)
# 0.005795514 0.005581524 0.005365759

# MSE back transformed (untransformed)
################################################################################
mse_untransformed_1val <- rowSums((imputed_widths_untransformed_1val -
                                       matrix(spruce_window[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                       length(impute_idx)
# 3728.496 3804.070 4202.224

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################
# Simulation study 4: Standard deviation is known
################################################################################

# We compare the AR(1), AR(2) and ARMA(1,1) with respect to how good they are
# at imputing missing values.

# One value to impute
#-------------------------------------------------------------------------------
imputed_widths_transformed_1val
convergence_mat # 0 means that optimization algorithm has converged


# Back transformation
# ------------------------------------------------------------------------------
t <- 1400:1800
coefs <- matrix(as.numeric(m$coef), nrow=2)
x <- cbind(rep(1, times=length(ts_res)), t)

log_val_mean <-x[impute_idx,]%*%coefs # alpha + beta*t

imputed_widths_untransformed_1val <- exp(imputed_widths_transformed_1val +
                                             matrix(log_val_mean, ncol=length(impute_idx),
                                                    nrow=3,byrow=T))*(rwl_sd_window_centered[impute_idx]+mean_sd)

# Transformed width, i.e before back transformation
# -------------------------------------------------
imputed_widths_transformed_1val_all <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_transformed_1val_all[,impute_idx] <- imputed_widths_transformed_1val
imputed_widths_transformed_1val_all[,impute_idx-1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_transformed_1val_all[,impute_idx+1] <- matrix(rep(ts_res,3), nrow=3, byrow=T)[,impute_idx+1]

# Width after back transformation
# -------------------------------
imputed_widths_1val <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_1val[,impute_idx] <- imputed_widths_untransformed_1val
imputed_widths_1val[,impute_idx-1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx-1]
imputed_widths_1val[,impute_idx+1] <- matrix(rep(spruce_window,3), nrow=3, byrow=T)[,impute_idx+1]

# True values
# ----------
true_values_transformed_1val <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val[-impute_idx] <- ts_res[-impute_idx]
true_values_untransformed_1val <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val[-impute_idx] <- spruce_window[-impute_idx]

# This contains the true values at the points where the imputatiions are done
true_values_transformed_1val_2 <- matrix(NA, nrow=1, ncol= length(ts_res))
true_values_transformed_1val_2[impute_idx] <- ts_res[impute_idx]
true_values_transformed_1val_2[impute_idx-1] <- ts_res[impute_idx-1]
true_values_transformed_1val_2[impute_idx+1] <- ts_res[impute_idx+1]

true_values_untransformed_1val_2 <- matrix(NA, nrow=1, ncol=length(ts_res))
true_values_untransformed_1val_2[impute_idx] <- spruce_window[impute_idx]
true_values_untransformed_1val_2[impute_idx-1] <- spruce_window[impute_idx-1]
true_values_untransformed_1val_2[impute_idx+1] <- spruce_window[impute_idx+1]

# Plots: transformed time series
plot(y=true_values_transformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_transformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_transformed_1val_all[1,], x=1:length(ts_res), col="green")
lines(imputed_widths_transformed_1val_all[2,], x=1:length(ts_res), col="blue")
lines(imputed_widths_transformed_1val_all[3,], x=1:length(ts_res), col="red")

# Plots: Untransformed (i.e. initial) time series
plot(y=true_values_untransformed_1val, x=1:length(ts_res), type="l", col="black")
lines(true_values_untransformed_1val_2[1,], x=1:length(ts_res), col="black", lty=3)
lines(imputed_widths_1val[1,], x=1:length(ts_res), col="green", lty=2)
lines(imputed_widths_1val[2,], x=1:length(ts_res), col="blue", lty=2)
lines(imputed_widths_1val[3,], x=1:length(ts_res), col="red", lty=2)

# MSE transformed
################################################################################
mse_transformed_1val <- rowSums((imputed_widths_transformed_1val -
                                     matrix(ts_res[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                     length(impute_idx)
# 0.005795514 0.005581524 0.005365759

# MSE back transformed (untransformed)
################################################################################
mse_untransformed_1val <- rowSums((imputed_widths_untransformed_1val -
                                       matrix(spruce_window[impute_idx], nrow=3, ncol=length(impute_idx), byrow=T))^2)/
                                       length(impute_idx)
# 4004.691 3903.181 4928.301

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


################################################################################
# Simulation study 4: Using na.kalman from the imputeTS package to impute the
# missing value
################################################################################

