source("make_series_stationary.R")

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

