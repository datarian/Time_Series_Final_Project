################################################################################

# Rk: the residuals are taken from the fourth approach

################################################################################


#library(MASS)
#mvrnorm(1, mu=c(0,0), Sigma=matrix(c(1,1,1,1),2,2))

rm(list = ls())

# Perhaps these files have to be run by hand
source("analysis/prepare_individual_series.R")
source("analysis//make_series_stationary.R")

# Literature
# ------------------------------------------------------------------------------
# https://sites.ualberta.ca/~sfossati/e509/files/other/dlm_ex.R
# https://robjhyndman.com/talks/ABS3.pdf


ts_res <- m$residuals
mean(ts_res) #4.697421e-18

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Imuputations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# We formulate an AR(1), AR(2) and an ARMA(1,1) as state space models. Imputations
# are done with the Kalman smoother. We also check the results with the Kalman
# smoother of th imputeTS package


# We also need to estimate rwl_sd_window
plot(y=rwl_sd_window, x=1:length(rwl_sd_window), type="l")
acf(rwl_sd_window)
pacf(rwl_sd_window) # AR(2)

mean_sd <- mean(rwl_sd_window)
rwl_sd_window_demeaned <- as.numeric(rwl_sd_window - mean_sd)


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
    sd <- x[2] #

    FF <- 1
    V <- 0
    GG <- phi
    W <- sd^2

    C0 <- sd^2/(1-phi^2)

    mAR1 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=1, C0=C0)
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
    V <- 1e-7
    W <- diag(c(sd^2, 0))

    mAR2 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=matrix(c(0,0), ncol=1), C0=diag(c(1, 1)))
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
# x_t = |phi 1|
#       |0   0|


buildARMA11 <- function(x){
    phi <- x[1]
    theta <- x[2]
    sd <- x[3]
    FF <- matrix(c(1,0), nrow=1)
    GG <- matrix(c(phi,0,1,0), 2,2)

    V <- 0
    K_W <-matrix(c(sd,1,1,theta*sd), 2,2)
    W <- K_W%*%t(K_W)

    mARMA11 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=matrix(c(0,0), ncol=1), C0=diag(c(1, 1)))
    return(mARMA11)
}


# We also need to estimate the variance. We model it as an AR(2) process

# Ex. start values
################################################################################
# Rk: ts_res has already mean 0

# AR(1)
#-------------------------------------------------------------------------------
AR1 <-arima(ts_res, order=c(1,0,0), include.mean = F)
phi_AR1_start <- as.numeric(AR1$coef) # parameter 1 to be optimized
#residuals_sd <- sqrt(400/401*var(ar1$residuals)) # parameter2 to be optimized
sigma_sd <- sqrt(as.numeric(AR1$sigma2))


# Simulation study
################################################################################

# We compare the AR(1), AR(2) and ARMA(1,1) with respect to how good they are
# at imputing missing values.

# One value to impute
#-------------------------------------------------------------------------------
set.seed(123)
impute_idx <- sample(1:401, 25)
impute_idx <- sort(impute_idx)


imputed_widths_untransformed_1val <- matrix(NA, nrow=3, ncol=length(impute_idx))
imputed_widths_transformed_1val <- imputed_widths_untransformed_1val
imputed_sd_1val <- rep(NA, length(impute_idx))

convergence_mat <- matrix(NA, nrow=4, ncol=length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- ts_res
    y_sd <- rwl_sd_window_demeaned
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
    fit <- dlmMLE(y, parm=initial_params, build=buildAR1, control = list(maxit = 1000))
    smoothed_state <- dlmSmooth(y, buildAR1(fit$par))$s[-1]
    imputed_widths_transformed_1val[1,j] <- smoothed_state[i] # F = 1
    convergence_mat[1,j] <- fit$convergence

    # AR(2)
    #...........................................................................
    AR2 <-arima(y, order=c(2,0,0), include.mean = F) # We assume that the series has still mean 0
    phi1 <- as.numeric(AR2$coef[1])
    phi2 <- as.numeric(AR2$coef[2])
    sigma <- sqrt(as.numeric(AR2$sigma2))

    initial_params <- c(phi1, phi2, sigma)
    fit <- dlmMLE(y, parm=initial_params, build=buildAR2, control = list(maxit = 1000))
    smoothed_states <- dlmSmooth(y, buildAR2(fit$par))$s
    smoothed_state1 <- smoothed_states[-1,1][i]
    smoothed_state2 <- smoothed_states[-1,2][i]
    imputed_widths_transformed_1val[2,j] <- matrix(c(1,0), nrow=1) %*%
        matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
    convergence_mat[2,j] <- fit$convergence

    # ARMA(1,1)
    #...........................................................................
    ARMA11 <-arima(y, order=c(1,0,1), include.mean = F) # We assume that the series has still mean 0
    phi <- as.numeric(ARMA11$coef[1])
    theta <- as.numeric(ARMA11$coef[2])
    sigma <- sqrt(as.numeric(ARMA11$sigma2))

    initial_params <- c(phi, theta, sigma)
    fit <- dlmMLE(y, parm=initial_params, build=buildARMA11, control = list(maxit = 1000))
    smoothed_states <- dlmSmooth(y, buildARMA11(fit$par))$s
    smoothed_state1 <- smoothed_states[-1,1][i]
    smoothed_state2 <- smoothed_states[-1,2][i]
    imputed_widths_transformed_1val[3,j] <- matrix(c(1,0), nrow=1) %*%
        matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
    convergence_mat[3,j] <- fit$convergence

    # Tree ring sd
    #---------------------------------------------------------------------------
    AR2 <-arima(y_sd, order=c(2,0,0), include.mean = F) # We assume that the series has still mean 0
    phi1 <- as.numeric(AR2$coef[1])
    phi2 <- as.numeric(AR2$coef[2])
    sigma <- sqrt(as.numeric(AR2$sigma2))

    initial_params <- c(phi1, phi2, sigma)
    fit <- dlmMLE(y, parm=initial_params, build=buildAR2, control = list(maxit = 1000))
    smoothed_states <- dlmSmooth(y, buildAR2(fit$par))$s
    smoothed_state1 <- smoothed_states[-1,1][i]
    smoothed_state2 <- smoothed_states[-1,2][i]
    imputed_sd_1val[j] <- matrix(c(1,0), nrow=1) %*%
        matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
    convergence_mat[4,j] <- fit$convergence

    j <- j + 1
}


imputed_widths_transformed_1val
imputed_sd_1val
convergence_mat


# Back transformation
# ------------------------------------------------------------------------------
t <- 1400:1800
#t2 <- t^2
coefs <- matrix(as.numeric(m$coef), nrow=2)
x <- cbind(rep(1, times=length(ts_res)), t)

log_val_mean <-x[impute_idx,]%*%coefs

imputed_widths_untransformed_1val <- exp(imputed_widths_transformed_1val +
                                         matrix(log_val_mean, ncol=length(impute_idx),
                                                nrow=3,byrow=T))*(imputed_sd_1val+mean_sd)

# Untransformed width, i.e before back transformation
imputed_widths_transformed_1val_all <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_transformed_1val_all[,impute_idx] <- imputed_widths_transformed_1val

# Width after back transformation
imputed_widths_1val <- matrix(NA, nrow=3, ncol=length(ts_res))
imputed_widths_1val[,impute_idx] <- imputed_widths_untransformed_1val


# Plots: transformed time series
plot(y=ts_res, x=1:length(ts_res), type="l")
points(imputed_widths_transformed_1val_all[1,], x=1:length(ts_res), col="green")
points(imputed_widths_transformed_1val_all[2,], x=1:length(ts_res), col="blue")
points(imputed_widths_transformed_1val_all[3,], x=1:length(ts_res), col="red")

# Plots: Untransformed (i.e. initial) time series
plot(y=spruce_window, x=1:length(ts_res), type="l")
points(imputed_widths_1val[1,], x=1:length(ts_res), col="green")
points(imputed_widths_1val[2,], x=1:length(ts_res), col="blue")
points(imputed_widths_1val[3,], x=1:length(ts_res), col="red")

# MSE untransformed
################################################################################
mse_transformed_1val <- rowSums((imputed_widths_transformed_1val -
                                     ts_res[impute_idx])^2)/length(impute_idx)
# 0.01528455 0.01473830 0.01522996

# MSE back transformed
################################################################################
mse_transformed_1val <- rowSums((imputed_widths_untransformed_1val -
                                     spruce_window[impute_idx])^2)/length(impute_idx)
# 996.3532 1083.8653  944.1404


