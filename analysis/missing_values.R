################################################################################

# Rk: the residuals are taken from the fourth approach

################################################################################

library(forecast)
library(matrixcalc)

# Goal and Approach
################################################################################
# Given the transformed times series of tree ring width data obtained by the
# fourth approach (see script 'make_series_stationary.R') we set 40 times randomly
# one observation as missing. Our goal is to 'predict' the missing value
# using all n-1 observation before and after the missing value. After having
# estimated the missing value we transform back the series to get the original
# time series. In order to do that we also have to estimate the process of the
# standard deviation (see 'make_series_stationary.R'). We calculate the
# mean-squared errors to compare the godness of fit. In particular, we model the
# transformed tree ring series either as AR(1), AR(2) or ARMA(1,1). The standard
# deviation is estimated as an AR(2) and as an ARMA(1,1) process. We also compare
# the results to the case when the standard deviation is known or linearly
# interpolated. All models are compared to the benchmark which is to use linear
# interpolation on the original time series.

# Literature:
# https://sites.ualberta.ca/~sfossati/e509/files/other/dlm_ex.R
# https://robjhyndman.com/talks/ABS3.pdf


# Missing values
################################################################################
# We divide the time series in 40 blocks of 9 data points. For our simulation
# study we will iterate through the 40 blocks and set randomly ONE observation of
# the actual block as missing. Hence, we will simulate 40 times the case, where
# ONE observation is missing and all the others are know. In order be able to
# plot the results in one graphics, the blocks are separated by one value. This
# ensures that two missing values never occur in sequence. In particular, the
# situation is like this:
# | 1, 2,..., 9| 10 | 11, 13, ..., 29 | 30 | 31, ..., 39| ... | 390 | 391, ..., 399| 400, 4001|
#     Block 1              Block 2                                       Block 40

# AR(1), AR(2) and ARMA(1,1) as state space models
################################################################################

# AR(1)
#-------------------------------------------------------------------------------
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
# Let x_t = [y_t, y_t-1]^T, and w_t = [e_t,0] and v_t = 0
# y_t = [1 0]x_t + v_t
# x_t = |phi_1 phi_2| x_{t-1} + w_t
#       |1        0 |


buildAR2 <- function(x){

    phi1 <- x[1]
    phi2 <- x[2]
    sd <- x[3]

    FF <- matrix(c(1,0), nrow=1)
    GG <- matrix(c(phi1, phi2, 1, 0), byrow=T, nrow=2)
    V <- 0
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

# K_W has the form matrix(c(1,theta,0,0), 2,2)*diag(sd^2, sd^2)
# Check for sd = 1
library(MASS)
theta <- -0.8
res <-mvrnorm(1, mu=c(0,0), Sigma=matrix(c(1,theta,theta,theta^2),2,2))
b <- res[1]
c <- res[2]
c/b

K <- matrix(c(1,theta,0,0), 2,2)
W <- K%*%t(K)
res <-mvrnorm(1, mu=c(0,0), Sigma=W)
b <- res[1]
c <- res[2]
c/b

svd(W)

buildARMA11 <- function(x){
    phi <- x[1]
    theta <- x[2]

    sd <- x[3]
    FF <- matrix(c(1,0), nrow=1)
    GG <- matrix(c(phi,0,1,0), 2,2)

    V <- 0

    K_W <- matrix(c(1,theta,0,0), 2,2) # K_W SHOULD BE LIKE THIS, BUT ...
    K_W <- K_W + matrix(c(1e-5,1e-5,1e-5,1e-5),2,2) # ... WE NEED TO ADD SOME NOISE ...
    W <- (K_W%*%t(K_W))%*%diag(c(sd^2, sd^2)) # ... OTHERWISE DLM WARNS THAT W IS ...
    # ... NOT A COVARIANCE MATRIX (IN PARTICULAR WHEN WE WANT TO MODEL THE SD TIME SERIES.

    mARMA11 <- dlm(FF=FF, GG=GG, V=V, W=W, m0=matrix(c(0,0), ncol=1), C0=diag(c(1e4, 1e4)))
    return(mARMA11)
}


# Transformed time series from the fourth aproach (see 'make_series_stationary.R')
################################################################################
ts_res <- m$residuals
mean(ts_res) #4.697421e-18
auto.arima(ts_res, d=0) # => ARIMA(1,0,1)


# Process of the standard deviation (see 'make_series_stationary.R')
################################################################################
mean_sd <- mean(rwl_sd_window)
rwl_sd_window_centered <- rwl_sd_window - mean_sd
auto.arima(rwl_sd_window_centered, d=0) # => ARIMA(1,0,1)

acf(rwl_sd_window_centered)
pacf(rwl_sd_window_centered) # => AR(2)


# Simulations
################################################################################

# Imputations on the transformed tree ring width time series
# ------------------------------------------------------------------------------
# Randomly chosen indices
set.seed(123)
impute_idx <- sample(1:9, 40, replace=T) + seq(0, 390, by=10)

# Matrices to store the imputed values and whether the optimization has converged.
# First row: AR(1); second row AR(2); third row: ARMA(1,1). 0 means 'converged'
# (see the dlm vignette).
imputed_widths_transformed <- matrix(NA, nrow=4, ncol=length(impute_idx))
convergence_mat_width <- matrix(NA, nrow=3, ncol=length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- ts_res
    y[i] <- NA

    # Imputations
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
        imputed_widths_transformed[1,j] <- smoothed_state[i] # F = 1
        convergence_mat_width[1,j] <- fit$convergence
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
        imputed_widths_transformed[2,j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat_width[2,j] <- fit$convergence
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
        convergence_mat_width[3,j] <- 1
    } else {
        smoothed_states <- dlmSmooth(y, buildARMA11(fit$par))$s
        smoothed_state1 <- smoothed_states[-1,1][i]
        smoothed_state2 <- smoothed_states[-1,2][i]
        imputed_widths_transformed[3,j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat_width[3,j] <- fit$convergence
    }

    # Linear interpolation
    # --------------------------------------------------------------------------
    imputed_widths_transformed[4,j] <- (y[i-1] + y[i+1])/2

    j <- j + 1
}


# Check of the results
# ------------------------------------------------------------------------------
imputed_widths_transformed
convergence_mat_width # 0 means that optimization algorithm has converged


# Imputations on the centered time series of the standard deviation
# ------------------------------------------------------------------------------
# Matrices to store the imputed values and whether the optimization has converged.
# First row: AR(2); second row: ARMA(1,1); third row: linear interpolation; fourth
# row: true standard deviation
# 0 means 'converged' (see the dlm vignette).
imputed_sd <- matrix(NA, nrow=4, ncol=length(impute_idx))
convergence_mat_sd <- matrix(NA, nrow=2, ncol=length(impute_idx))


j <- 1
for (i in impute_idx){

    print(j)

    y_sd <- rwl_sd_window_centered
    y_sd[i] <- NA

    # Imputations
    ############################################################################

    # Tree ring sd as AR(2)
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
        imputed_sd[1,j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat_sd[1,j] <- fit$convergence
    }


    # Tree ring sd as ARMA(1,1)
    #---------------------------------------------------------------------------
    ARMA11 <-arima(y_sd, order=c(1,0,1), include.mean = F) # We assume that the series has still mean 0
    phi <- as.numeric(ARMA11$coef[1])
    theta <- as.numeric(ARMA11$coef[2])
    sigma <- sqrt(as.numeric(ARMA11$sigma2))

    initial_params <- c(phi, theta, sigma)
    fit <- try(dlmMLE(y_sd, parm=initial_params, build=buildARMA11, control = list(maxit = 1000), method="BFGS"))
    if(class(fit)=="try-error"){
        convergence_mat[2,j] <- 1

    } else {
        smoothed_states <- dlmSmooth(y_sd, buildARMA11(fit$par))$s
        smoothed_state1 <- smoothed_states[-1,1][i]

        smoothed_state2 <- smoothed_states[-1,2][i]
        imputed_sd[2,j] <- matrix(c(1,0), nrow=1) %*%
            matrix(c(smoothed_state1, smoothed_state2), ncol=1) # F = [1 0]
        convergence_mat_sd[2,j] <-fit$convergence
    }

    # Tree ring sd linear interpolation
    #---------------------------------------------------------------------------
    imputed_sd[3,j] <- (y_sd[i-1] + y_sd[i+1])/2


    # True tree ring sd
    #---------------------------------------------------------------------------
    imputed_sd[4,j] <- rwl_sd_window_centered[i]

    j <- j + 1
}

imputed_sd
convergence_mat_sd # all entries are 0


# Coefficients of linear model for back-transformation
################################################################################
# Matrix to store the coefficients of the linear model. First row: intercept;
# Second row: slope
coefficients_lin_model_mat <- matrix(NA, nrow=2, ncol=length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- spruce_window
    y[i] <- NA

    y_sd <- rwl_sd_window
    y_sd[i] <- NA

    y_stab <- y/y_sd
    y_stab_log <- log(y_stab)

    # Linear model
    # --------------------------------------------------------------------------
    m <- lm(y_stab_log ~ t)
    coefficients_lin_model_mat[1,j] <- m$coefficients[1]
    coefficients_lin_model_mat[2,j] <- m$coefficients[2]

    j <- j + 1
}

coefficients_lin_model_mat

# Back-transformation
################################################################################
t <- 1400:1800
X <- cbind(rep(1, times=length(ts_res)), t)

log_val_mean_imputed_values <- X[impute_idx,]%*%coefficients_lin_model_mat
log_val_mean_imputed_values <- diag(log_val_mean_imputed_values)

# Add back the linear regression line and exponentiate
# ------------------------------------------------------------------------------
y_transformed_mat <- exp(matrix(log_val_mean_imputed_values, nrow=nrow(imputed_widths_transformed)
                                , ncol=length(impute_idx), byrow=T) + imputed_widths_transformed)


# Case 1: sd as AR(2)
# ------------------------------------------------------------------------------
y_imputed_sd_ar2_mat <- y_transformed_mat*matrix(imputed_sd[1,]+mean_sd, nrow=nrow(imputed_widths_transformed),
                                                 ncol=length(impute_idx), byrow=T)

# Case 2: sd as ARMA(1,1)
# ------------------------------------------------------------------------------
y_imputed_sd_arma11_mat <- y_transformed_mat*matrix(imputed_sd[2,]+mean_sd, nrow=nrow(imputed_widths_transformed),
                                                    ncol=length(impute_idx), byrow=T)

# Case 3: sd linearlly interpolated
# ------------------------------------------------------------------------------
y_imputed_sd_lin_int_mat <- y_transformed_mat*matrix(imputed_sd[3,]+mean_sd, nrow=nrow(imputed_widths_transformed),
                                                     ncol=length(impute_idx), byrow=T)

# Case 4: true sd
# ------------------------------------------------------------------------------
y_imputed_sd_true_mat <- y_transformed_mat*matrix(imputed_sd[4,]+mean_sd, nrow=nrow(imputed_widths_transformed),
                                                  ncol=length(impute_idx), byrow=T)

# Linear Interpolatin of the untransformed time series
################################################################################
# One of the simplest method is to linearlly interpolated the original time series
y_original_tmp1 <- spruce_window[impute_idx-1]
y_original_tmp2 <- spruce_window[impute_idx+1]

y_original_lin_int <- (y_original_tmp1 + y_original_tmp2)/2

# Model comparison
################################################################################

# MSE transformed time series
# ------------------------------------------------------------------------------
mse_transformed <- rowSums((imputed_widths_transformed -
                                matrix(ts_res[impute_idx], nrow=nrow(imputed_widths_transformed),
                                ncol=length(impute_idx), byrow=T))^2)/ length(impute_idx)

# 0.005795514 0.005581524 0.005365759 0.005660842 # ARIMA(1,1) is the best model, followed by AR(2) and ...
                                                  # linear interpolation. AR(1) is worst.

# MSE back-transformed time series
# ------------------------------------------------------------------------------
# sd as AR(2)
mse_back_transformed_sd_ar2 <- rowSums((y_imputed_sd_ar2_mat - matrix(spruce_window[impute_idx],
                                    nrow=nrow(imputed_widths_transformed),
                                    ncol=length(impute_idx), byrow=T))^2)/length(impute_idx)

# 510.1684 492.8501 484.3971 527.1005 # Linear Interpolation of the tree ring width is worst

# sd as ARMA(1,1)
mse_back_transformed_sd_arma11 <- rowSums((y_imputed_sd_arma11_mat - matrix(spruce_window[impute_idx],
                                    nrow=nrow(imputed_widths_transformed),
                                    ncol=length(impute_idx), byrow=T))^2)/length(impute_idx)

# 286.4674 281.7908 274.1983 306.1778 # Linear Interpolation of the tree ring width is worst

# sd linearly interpolated
mse_back_transformed_sd_lin_int <- rowSums((y_imputed_sd_lin_int_mat - matrix(spruce_window[impute_idx],
                                    nrow=nrow(imputed_widths_transformed),
                                    ncol=length(impute_idx), byrow=T))^2)/length(impute_idx)

# 273.5193 278.2707 275.7508 285.8944 # Linear Interpolation of the tree ring width is worst

# true sd
mse_back_transformed_sd_true <- rowSums((y_imputed_sd_true_mat - matrix(spruce_window[impute_idx],
                                    nrow=nrow(imputed_widths_transformed),
                                    ncol=length(impute_idx), byrow=T))^2)/length(impute_idx)

# 175.8117 170.2908 163.2709 169.1063 # Linear interpolation is second

# MSE of the time series resulting from linearlly interpolate the original time
# series
# ------------------------------------------------------------------------------
mes_original_lin_int <- mean((y_original_lin_int-spruce_window[impute_idx])^2)
# 302.3888 This is slightly better than assuming ARMA(1,1) for width and sd, but worse than ...
#          assuming ARMA(1,1) for the width and linear interpolation for the sd

# Summary
# ------------------------------------------------------------------------------
# Assuming an ARMA(1,1) process as the process of the transformed tree ring width
# gives the best results when back-transformed no matter how the centered standard
# deviation is calculated. The ARMA(1,1) process is also best compared to the
# transformed tree ring width time series.
# Imputations knowing the true variance is best as expected. However, simple linear
# interpolation outperforms the case where we assumed an AR(2) or an ARMA(1,1),
# respectively as the process of the centered sd.
# Assuming ARMA(1,1) for the process of the transformed tree ring width and ARMA(1,1)
# or linear interpolation for the process of the centered standard deviation
# outperforms all the other approaches.


# Check with the imputeTS package
################################################################################

imputed_widths_transformed_imputeTS <- numeric(length(impute_idx))
imputed_sd_imputeTS <- numeric(length(impute_idx))

j <- 1
for (i in impute_idx){

    print(j)

    y <- ts_res
    y[i] <- NA

    y_sd <- rwl_sd_window_centered
    y_sd[i] <- NA

    # Imputations
    ############################################################################

    tmp <- na.kalman(y, model="auto.arima")
    imputed_widths_transformed_imputeTS[j] <- tmp[i]

    tmp <- na.kalman(y_sd, model="auto.arima")
    imputed_sd_imputeTS[j] <- tmp[i]

    j <- j + 1
}

# MSE own implementation vs. imputeTS package
# ------------------------------------------------------------------------------
mse_own_package_transformed_width <- mean((imputed_widths_transformed_imputeTS -
                                              imputed_widths_transformed[3,])^2)
# 8.814445e-05

mse_own_package_sd <- mean((imputed_sd_imputeTS - imputed_sd[2,])^2)
# 2.660028

# Plots own implementation vs. imputeTS package
# ------------------------------------------------------------------------------
# Widths
plot(y=imputed_widths_transformed[3,], x=1:length(impute_idx), type="p")
points(y=imputed_widths_transformed_imputeTS, x=1:length(impute_idx), col="blue")

# Centered sd
plot(y=imputed_sd[2,], x=1:length(impute_idx), type="p")
points(y=imputed_sd_imputeTS, x=1:length(impute_idx), col="blue")

# As ggplots
# Widths
#.......
df_widths <- data.frame(Transformed_Width=c(imputed_widths_transformed[3,], imputed_widths_transformed_imputeTS),
                        Index = c(impute_idx, impute_idx),
                        Method = c(rep(1, length(impute_idx)), rep(2, length(impute_idx))))
df_widths$Method <- factor(df_widths$Method, labels=c("Own implementation", "imputeTS"))


pplot_width_comparison <- ggplot(df_widths, aes(y=Transformed_Width, x=Index, group=Method))
pplot_width_comparison + geom_point(aes(col=Method), alpha = 0.8) +
    scale_color_manual(values=c("blue", "green"))+
    labs(title="Imputed widths (transformed series): own implementation vs. imputeTS",
         y="Width") +
    theme(legend.position="bottom",
          plot.title = element_text(face="bold", hjust=0.5))

# Centered sd
#............
df_sd <- data.frame(Centered_SD=c(imputed_sd[2,], imputed_sd_imputeTS),
                        Index = c(impute_idx, impute_idx),
                        Method = c(rep(1, length(impute_idx)), rep(2, length(impute_idx))))
df_sd$Method <- factor(df_widths$Method, labels=c("Own implementation", "imputeTS"))


pplot_width_comparison <- ggplot(df_sd, aes(y=Centered_SD, x=Index, group=Method))
pplot_width_comparison + geom_point(aes(col=Method), alpha = 0.8) +
    scale_color_manual(values=c("blue", "green"))+
    labs(title="Imputed SD (centered series): own implementation vs. imputeTS",
         y="Centered SD") +
    theme(legend.position="bottom",
          plot.title = element_text(face="bold", hjust=0.5))


# Is linear interpolation for the sd still as good as assuming an ARMA(1,1)?
# ------------------------------------------------------------------------------
# Here we check if it is still equivalent in terms of MSE to linearly interpolate
# the sd. The alternative is to assume an ARMA(1,1) process. The process of the
# transformed tree ring width is ARMA(1,1).

# Back-transformation
# ..............................................................................
y_transformed_imputeTS <- exp(matrix(log_val_mean_imputed_values, nrow=1,
                                ncol=length(impute_idx), byrow=T) + imputed_widths_transformed_imputeTS)

y_back_transformed_sd_arma11_imputeTS <- y_transformed_imputeTS*(imputed_sd_imputeTS+mean_sd)
y_back_transformed_sd_lin_int_imputeTS <- y_transformed_imputeTS*(imputed_sd[3,]+mean_sd)

# MSE back-transformed vs. linearlly interpolated sd vs. linearlly interpolated orginal values
# ..............................................................................
mse_sd_arma11_imputeTS <- mean((spruce_window[impute_idx]-y_back_transformed_sd_arma11_imputeTS)^2)
# 256.1351 # This value is sligthly better than with our implementation

mse_sd_lin_int_imputeTS <- mean((spruce_window[impute_idx]-y_back_transformed_sd_lin_int_imputeTS)^2)
# 275.6912 # Same value as with our implementation


# Plot of the time series with the imputed values
################################################################################
# Here we plot the resulting time series using our implementation and  the imputeTS
# package. We also plot the results that we achieve by simply linearly interpolate
# the orginal time series.

# Own implementation
# ------------------------------------------------------------------------------
own_implementation <- rep(NA, times=length(spruce_window))
own_implementation[impute_idx] <- y_imputed_sd_arma11_mat[3,]
if (impute_idx[1] > 1){
    own_implementation[impute_idx-1] <- spruce_window[impute_idx-1]
} else {
    own_implementation[impute_idx[-1]-1] <- spruce_window[impute_idx[-1]-1] # first index could be 1
}
own_implementation[impute_idx+1] <- spruce_window[impute_idx+1]

# ImputeTS
# ------------------------------------------------------------------------------
imputeTS <- rep(NA, times=length(spruce_window))
imputeTS[impute_idx] <- y_back_transformed_sd_arma11_imputeTS
if (impute_idx[1] > 1){
    imputeTS[impute_idx-1] <- spruce_window[impute_idx-1]
} else {
    imputeTS[impute_idx[-1]-1] <- spruce_window[impute_idx[-1]-1] # first index could be 1
}
imputeTS[impute_idx+1] <- spruce_window[impute_idx+1]

# Linear interpolation of the original time series
# ------------------------------------------------------------------------------
lin_int <- rep(NA, times=length(spruce_window))
lin_int[impute_idx] <- as.numeric(y_original_lin_int)
if (impute_idx[1] > 1){
    lin_int[impute_idx-1] <- spruce_window[impute_idx-1]
} else {
    lin_int[impute_idx[-1]-1] <- spruce_window[impute_idx[-1]-1] # first index could be 1
}
lin_int[impute_idx+1] <- spruce_window[impute_idx+1]

# True values with NA where we impute
# ------------------------------------------------------------------------------
true_without_imputed <- spruce_window
true_without_imputed[impute_idx] <- NA

# True values
# ------------------------------------------------------------------------------
true_imputed <- rep(NA, times=length(spruce_window))
true_imputed[impute_idx] <- as.numeric(spruce_window[impute_idx])
if (impute_idx[1] > 1){
    true_imputed[impute_idx-1] <- spruce_window[impute_idx-1]
} else {
    true_imputed[impute_idx[-1]-1] <- spruce_window[impute_idx[-1]-1] # first index could be 1
}
true_imputed[impute_idx+1] <- spruce_window[impute_idx+1]


# Data frame for the ggplot
# ------------------------------------------------------------------------------
data <- data.frame(Values=c(own_implementation, imputeTS, lin_int, true_imputed, true_without_imputed),
                            Time = rep(1400:1800, 5),
                            Method = c(rep(1, length(spruce_window)), rep(2, length(spruce_window)),
                                           rep(3, length(spruce_window)), rep(4, length(spruce_window)),
                                       rep(5, length(spruce_window))))

data$Method <- factor(data$Method, labels=c("Own implementation", "imputeTS",
                                    "Linear interpolation", "True imputation", "Original series"))

data_NA_removed <- data[!is.na(data$Values),]
tail(data_NA_removed)

pplot_comparison <- ggplot(data, aes(y=Values, x=Time, group=Method))
pplot_comparison <- pplot_comparison + geom_line(aes(col=Method), alpha = 1) +
    scale_color_manual(values=c("blue", "green", "purple", "red", "black")) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "solid")) +
    theme(legend.position="bottom",
          plot.title = element_text(face="bold", hjust=0.5)) +
    geom_vline(xintercept = t[impute_idx], col="grey")

# Summary
################################################################################
# The results with the imputeTS package are slighty better. In particlar, the
# imputaions of the centered sd in case of assuming an ARMA(1,1) process works
# slightly better than with our implementation. A reason could be that we have
# to add some noise to the covariance matrix when implementing the ARMA(1,1) process
# as state space model so that the dlm package is able to do the cholesky decomposition.
# The imputeTS package seems to handle this better.



