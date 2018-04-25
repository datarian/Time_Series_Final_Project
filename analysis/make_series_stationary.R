# getting window of the two series from 1332 to 1800:
spruce_window <- window(spruce_sup_900_ts,end=c(1800))

log_spruce_window <- log(spruce_window)
t <- seq(1,length(spruce_window))
t2 <- t^2

d <- data.frame(cbind(log_spruce_window,t,t2))



m2 <- lm(log_spruce_window ~ t+t2,data=d)
summary(m2)

# Transformation of residuals
i <- (m2$residuals/m2$fitted.values)
