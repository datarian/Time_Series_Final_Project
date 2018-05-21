
# Finding candidate time series models and assessing them.



Xt <- order1_model$residuals # the stationary series created in make_series_stationary.R
Xt.bar <- mean(Xt)

# Find best model automatically:
Xt.autofit <- autofit(x=Xt,p=1:3,q=0:3) # ARMA(1,1) was found.


#*******************************************************************************
# Fit 3 models: ARMA(1,1), AR(1), AR(2) and extract parameters for comparison
# Decision basis: literature hints, ACF/PACF interpretaion


# Extract parameters and model residuals

# Arma(1,1)
Xt.arma11 <- arima(x=Xt,order=c(1,0,1),include.mean = F)
arma11.phi <- Xt.arma11$model$phi
arma11.phi.causal <- abs( polyroot(c(1,-arma11.phi)))
arma11.theta <- Xt.arma11$model$theta
arma11.sigma2 <- Xt.arma11$sigma2
arma11_resids <- cbind(year=t,residuals=Xt.arma11$residuals/sd(Xt.arma11$residuals),model=rep("ARMA(1,1)",times=length(t)))

# AR(1)
Xt.ar1 <- arima(x=Xt,order=c(1,0,0),include.mean = F)
ar1.phi <- Xt.ar1$model$phi
ar1.phi.causal <- abs( polyroot(c(1,-ar1.phi)))
ar1.sigma2 <- Xt.ar1$sigma2
ar1_resids <- cbind(year=t,residuals=Xt.ar1$residuals/sd(Xt.ar1$residuals),model=rep("AR(1)",times=length(t)))

# AR(2)
Xt.ar2 <- arima(x=Xt,order=c(2,0,0),include.mean = F)
ar2.phi <- Xt.ar2$model$phi
ar2.phi.causal <- abs( polyroot(c(1,-ar2.phi)))
ar2.sigma2 <- Xt.ar2$sigma2
ar2_resids <- cbind(year=t,residuals=Xt.ar2$residuals/sd(Xt.ar2$residuals),rep("AR(2)",times=length(t)))


#*******************************************************************************
# Prepare data frames for report tables


# Compare AIC's of the 3 models
model_aic_compared <- data.frame(Model=c("ARMA(1,1)","AR(1)","AR(2)"),
                                 AIC=c(Xt.arma11$aic,Xt.ar1$aic,Xt.ar2$aic))


# Build a comparison dataframe for parameters
paramnames <- c("$\\phi$","$\\theta$","$\\sigma^2_{ARMA(1,1)}$",
                "$\\phi$","$\\sigma^2_{AR(1)}$",
                "$\\phi_1$","$\\phi_2$","$\\sigma^2_{AR(2)}$")

parameters <- round(c(arma11.phi,arma11.theta,arma11.sigma2,
                      ar1.phi,ar1.sigma2,
                      ar2.phi,ar2.sigma2),digits=3)

se <- c(round(sqrt(diag(Xt.arma11$var.coef)),digits=3),"-",
        round(sqrt(diag(Xt.ar1$var.coef)),digits=3),"-",
        round(sqrt(diag(Xt.ar2$var.coef)),digits=3),"-")

modelnames <- c(rep("ARMA(1,1)",times=3),
                rep("AR(1)",times=2),
                rep("AR(2)",times=3))

modelComparisonTable <- data.frame(modelnames,paramnames,parameters,se)
colnames(modelComparisonTable) <- c("Model","Parameter name",
                                    "Parameter value", "standard error")


#*******************************************************************************
# Manually build up the tsdiag() plots for all three models combined

# ACF

acf_arma11 <- acf(Xt.arma11$residuals, type="correlation",plot=F)
acf_arma11 <- data.frame(lag=0:26,
                         acf=acf_arma11$acf,
                         model=rep("ARMA(1,1)",times=27))
acf_ar1 <- acf(Xt.ar1$residuals, type="correlation",plot=F)
acf_ar1 <- data.frame(lag=0:26,
                      acf=acf_ar1$acf,
                      model=rep("AR(1)",times=27))
acf_ar2 <- acf(Xt.ar2$residuals, type="correlation",plot=F)
acf_ar2 <- data.frame(lag=0:26,
                      acf=acf_ar2$acf,
                      model=rep("AR(2)",times=27))
acf_comparison_df <- rbind(
    acf_arma11,
    acf_ar1,
    acf_ar2
)
acf_comparison_df$lag <- as.factor(acf_comparison_df$lag)
acf_comparison_df$model <- factor(acf_comparison_df$model,
                                  levels = c("AR(1)","AR(2)","ARMA(1,1)"))

acf_comparison_plot <- ggplot(acf_comparison_df,aes(x=lag,y=acf,fill=model)) +
    geom_col(position="dodge", width = 0.7) +
    geom_hline(yintercept=0,size=0.3) +
    geom_hline(yintercept=c(-1.96/sqrt(400),1.96/sqrt(400)),colour="darkgrey",size=0.3) +
    ylab("ACF")


# Standardized residuals (no need to manually standardize, they are already scaled)

residuals_df <- as.data.frame(rbind(arma11_resids,ar1_resids,ar2_resids))
residuals_df$year <- as.numeric(levels(residuals_df$year)[residuals_df$year])
residuals_df$residuals <- as.numeric(levels(residuals_df$residuals)[residuals_df$residuals])
residuals_df$model <- factor(residuals_df$model, levels = c("AR(1)","AR(2)","ARMA(1,1)"))

residuals_plot <- ggplot(residuals_df, aes(x=year,y=residuals,colour=model)) +
    geom_hline(yintercept=0,size=0.3) +
    geom_line(size=0.3)


# Ljung-Box test

boxtest_df <- data.frame()

for(l in 1:10){
    if(l>1){
        p_box_ar1 <- Box.test(Xt.ar1$residuals,
                              lag=l,type="Ljung-Box",fitdf = 2)$p.value
        box_ar1 <- data.frame(lag=l, pvalue=p_box_ar1, model="AR(1)")

        boxtest_df <- rbind(boxtest_df, box_ar1)
    }
    if(l>2){
        p_box_arma11 <- Box.test(Xt.arma11$residuals,
                                 lag=l,type="Ljung-Box",fitdf = 2)$p.value
        box_arma11 <- data.frame(lag=l, pvalue=p_box_arma11, model="ARMA(1,1)")

        p_box_ar2 <- Box.test(Xt.ar2$residuals,
                              lag=l,type="Ljung-Box",fitdf = 2)$p.value
        box_ar2 <- data.frame(lag=l, pvalue=p_box_ar2, model="AR(2)")
        boxtest_df <- rbind(boxtest_df, box_arma11, box_ar2)
    }
}
colnames(boxtest_df) <- c("lag", "pvalue", "model")
boxtest_df$model <- factor(boxtest_df$model,
                           levels = c("AR(1)","AR(2)","ARMA(1,1)"))

boxtest_plot <- ggplot(boxtest_df, aes(x=lag,y=pvalue,colour=model)) +
    geom_point() +
    geom_hline(yintercept=0.05, colour="darkgrey",size=0.3) +
    ylab("p value")
