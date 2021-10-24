library(nlme)
library(car)

data <- read.delim("./data_102018.txt")

########################
### ALL CRASH CONTROLLED
########################

all_ols <- lm(tot_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m9, data=data)
summary(all_ols)

# Durbin-watson test, 12 time periods
dwt(all_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serially correlated errors
plot(data$time[1:192],
     residuals(all_ols)[1:192],
     type='o',
     pch=16,
     xlab='Time',
     ylab='OLS Residuals',
     col="red")
abline(h=0,lty=2)

# Plot ACF and PACF
# Set plotting to two records on one page
par(mfrow=c(1,2))

# Produce plots
acf(residuals(all_ols))
acf(residuals(all_ols),type='partial')
# Note decay in ACF, significant spike at 1 in PACF, model p=1

########################
# Run the final model
########################

# Fit the GLS regression model
all_p1 <- gls(tot_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m9,
              data=data,
              correlation=corARMA(p=1,form=~time|tsp),
              method="ML")
summary(all_p1)
confint(all_p1)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
all_p2 <- update(all_p1,correlation=corARMA(p=2,form=~time|tsp))
anova(all_p1,all_p2)

all_p1q1 <- update(all_p1,correlation=corARMA(p=1,q=1,form=~time|tsp))
anova(all_p1,all_p1q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(all_p1))


########################
# Plot results
#########################

# First plot the raw data points for the tsp corridors
par(mfrow=c(1,1))
plot(data$time[1:192],data$tot_l[1:192],
     ylim=c(0,250),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add in the points for the control
points(data$time[193:384],data$tot_l[193:384],
       col="lightblue",
       pch=20)

# Add x-axis year labels
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Label the deployment of tsp
abline(v=89.5,lty=2)

# Calculate the seasonality offsets
offset <- mean(data$m9) * all_p1$coef[9]

# Plot the first line segment of the treated group
segments(1,
         all_p1$coef[1] + all_p1$coef[3] + all_p1$coef[2] 
         + all_p1$coef[4] + offset,
         89,
         all_p1$coef[1] + all_p1$coef[3] 
         + (all_p1$coef[2] + all_p1$coef[4])*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment of the treated group
segments(90,
         all_p1$coef[1] + all_p1$coef[3] 
         + (all_p1$coef[2] + all_p1$coef[4])*90 
         + all_p1$coef[5] + all_p1$coef[6] 
         + all_p1$coef[7] + all_p1$coef[8] + offset,
         192,
         all_p1$coef[1] + all_p1$coef[3] 
         + (all_p1$coef[2] + all_p1$coef[4])*192
         + all_p1$coef[5] + all_p1$coef[7]
         + (all_p1$coef[6] + all_p1$coef[8])*103 + offset,
         lty=1, lwd=2, col='red')

# Add the counterfactual for the treated group
segments(90, all_p1$coef[1] + all_p1$coef[2]*90 
         + all_p1$coef[3] + all_p1$coef[4]*90 
         + all_p1$coef[5] + all_p1$coef[6] + offset,
         192, all_p1$coef[1] + all_p1$coef[2]*192 
         + all_p1$coef[3] + all_p1$coef[4]*90 
         + all_p1$coef[5] + all_p1$coef[6]*103 + offset,
         lty=2,col='red',lwd=2)

# Plot the first line segment of the control group
segments(1,
         all_p1$coef[1] + all_p1$coef[2] + offset,
         89,
         all_p1$coef[1] + all_p1$coef[2]*89 + offset,
         lty=1, lwd=2, col='blue')

# Plot the second line segment of the control group
segments(90,
         all_p1$coef[1] + all_p1$coef[2]*90 + 
           all_p1$coef[5] + all_p1$coef[6] + offset,
         192,
         all_p1$coef[1] + all_p1$coef[2]*192 + 
           all_p1$coef[5] + all_p1$coef[6]*103 + offset,
         lty=1, lwd=2, col='blue')

# Add the counterfactual for the control group
segments(1, all_p1$coef[1]+all_p1$coef[2] + offset,
         192,all_p1$coef[1]+all_p1$coef[2]*192 + offset,
         lty=2,col='blue',lwd=2)

# Add in a legend
legend("topright", inset=.02, legend=c("tsp","control"), col=c("red","blue"),pch=20)

##############################################
# Predict absolute and relative changes
##############################################

# Predicted crashes after tsp implementation

for (i in 1:103){
  pred <- fitted(all_p1)[89+i]
  print(unname(pred))
}

# Estimate the counterfactual for the same period

for (i in 1:103){
  cfac <- all_p1$coef[1] + all_p1$coef[2]*(89+i) +
          all_p1$coef[3] + all_p1$coef[4]*90 +
          all_p1$coef[5] + all_p1$coef[6]*i
  print(unname(cfac))
}

#END


########################
### PDO CRASH CONTROLLED
########################

pdo_ols <- lm(pdo_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m3+m6+m7+m9, data=data)
summary(pdo_ols)

# Durbin-watson test, 12 time periods
dwt(pdo_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for seripdoy correlated errors
plot(data$time[1:192],
     residuals(pdo_ols)[1:192],
     type='o',
     pch=16,
     xlab='Time',
     ylab='OLS Residuals',
     col="red")
abline(h=0,lty=2)

# Plot ACF and PACF
# Set plotting to two records on one page
par(mfrow=c(1,2))

# Produce plots
acf(residuals(pdo_ols))
acf(residuals(pdo_ols),type='partial')
# Note decay in ACF, significant spike at 1 in PACF, model p=1

########################
# Run the final model
########################

# Fit the GLS regression model
pdo_p1 <- gls(pdo_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m3+m6+m7+m9,
              data=data,
              correlation=corARMA(p=1,form=~time|tsp),
              method="ML")
summary(pdo_p1)
confint(pdo_p1)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
pdo_p2 <- update(pdo_p1,correlation=corARMA(p=2,form=~time|tsp))
anova(pdo_p1,pdo_p2)

pdo_p1q1 <- update(pdo_p1,correlation=corARMA(p=1,q=1,form=~time|tsp))
anova(pdo_p1,pdo_p1q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(pdo_p1))


########################
# Plot results
########################

# First plot the raw data points for the tsp corridors
par(mfrow=c(1,1))
plot(data$time[1:192],data$pdo_l[1:192],
     ylim=c(0,150),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add in the points for the control
points(data$time[193:384],data$pdo_l[193:384],
       col="lightblue",
       pch=20)

# Add x-axis year labels
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Label the deployment of tsp
abline(v=89.5,lty=2)

# Calculate the seasonality offsets
offset <- mean(data$m3) * pdo_p1$coef[9]
          +mean(data$m6) * pdo_p1$coef[10]
          +mean(data$m7) * pdo_p1$coef[11]
          +mean(data$m9) * pdo_p1$coef[12]

# Plot the first line segment of the treated group
segments(1,
         pdo_p1$coef[1] + pdo_p1$coef[3] + pdo_p1$coef[2] 
         + pdo_p1$coef[4] + offset,
         89,
         pdo_p1$coef[1] + pdo_p1$coef[3] 
         + (pdo_p1$coef[2] + pdo_p1$coef[4])*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment of the treated group
segments(90,
         pdo_p1$coef[1] + pdo_p1$coef[3] 
         + (pdo_p1$coef[2] + pdo_p1$coef[4])*90 
         + pdo_p1$coef[5] + pdo_p1$coef[6] 
         + pdo_p1$coef[7] + pdo_p1$coef[8] + offset,
         192,
         pdo_p1$coef[1] + pdo_p1$coef[3] 
         + (pdo_p1$coef[2] + pdo_p1$coef[4])*192
         + pdo_p1$coef[5] + pdo_p1$coef[7]
         + (pdo_p1$coef[6] + pdo_p1$coef[8])*103 + offset,
         lty=1, lwd=2, col='red')

# Add the counterfactual for the treated group
segments(90, pdo_p1$coef[1] + pdo_p1$coef[2]*90 
         +pdo_p1$coef[3] + pdo_p1$coef[4]*90
         + pdo_p1$coef[5] + pdo_p1$coef[6] + offset,
         192, pdo_p1$coef[1] + pdo_p1$coef[2]*192 
         + pdo_p1$coef[3] + pdo_p1$coef[4]*90 
         + pdo_p1$coef[5] + pdo_p1$coef[6]*103 + offset,
         lty=2,col='red',lwd=2)

# Plot the first line segment of the control group
segments(1,
         pdo_p1$coef[1] + pdo_p1$coef[2] + offset,
         89,
         pdo_p1$coef[1] + pdo_p1$coef[2]*89 + offset,
         lty=1, lwd=2, col='blue')

# Plot the second line segment of the control group
segments(90,
         pdo_p1$coef[1] + pdo_p1$coef[2]*90 + 
           pdo_p1$coef[5] + pdo_p1$coef[6] + offset,
         192,
         pdo_p1$coef[1] + pdo_p1$coef[2]*192 + 
           pdo_p1$coef[5] + pdo_p1$coef[6]*103 + offset,
         lty=1, lwd=2, col='blue')

# Add the counterfactual for the control group
segments(1, pdo_p1$coef[1] + pdo_p1$coef[2] + offset,
         192,pdo_p1$coef[1] + pdo_p1$coef[2]*192 + offset,
         lty=2,col='blue',lwd=2)

# Add in a legend
legend("topright", inset=.02, legend=c("tsp","control"), col=c("red","blue"),pch=20)


##############################################
# Predict absolute and relative changes
##############################################

# Predicted crashes after tsp implementation
#pred <- 0
#for (i in 1:103){
#  pred <- pred + fitted(pdo_p1)[89+i]
#}

# Estimate the counterfactual for the same period
#cfac <- 0
#for (i in 1:103){
#  cfac <- (cfac + pdo_p1$coef[1] + pdo_p1$coef[2]*(89+i) +
#             pdo_p1$coef[3] + pdo_p1$coef[4]*(89+i) +
#             pdo_p1$coef[5] + pdo_p1$coef[6]*i)
#}

#pred
#cfac

# Absolute change after tsp implementation
#pred - cfac
# Relative change after tsp implementation
#(pred - cfac) / cfac

## evaluation::

# Predicted crashes after tsp implementation

for (i in 1:103){
  pred <- fitted(pdo_p1)[89+i]
  print(unname(pred))
}

# Estimate the counterfactual for the same period

for (i in 1:103){
  cfac <- pdo_p1$coef[1] + pdo_p1$coef[2]*(89+i) +
          pdo_p1$coef[3] + pdo_p1$coef[4]*90 +
          pdo_p1$coef[5] + pdo_p1$coef[6]*i
  print(unname(cfac))
}

#END



########################
### FI CRASH CONTROLLED
########################

fi_ols <- lm(fi_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m5+m8, data=data)
summary(fi_ols)

# Durbin-watson test, 12 time periods
dwt(fi_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serifiy correlated errors
plot(data$time[1:192],
     residuals(fi_ols)[1:192],
     type='o',
     pch=16,
     xlab='Time',
     ylab='OLS Residuals',
     col="red")
abline(h=0,lty=2)

# Plot ACF and PACF
# Set plotting to two records on one page
par(mfrow=c(1,2))

# Produce plots
acf(residuals(fi_ols))
acf(residuals(fi_ols),type='partial')
# Note decay in ACF, significant spike at 5 in PACF, model p=5

########################
# Run the final model
########################

# Fit the GLS regression model
fi_p5 <- gls(fi_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m5+m8,
              data=data,
              correlation=corARMA(p=5,form=~time|tsp),
              method="ML")
summary(fi_p5)
confint(fi_p5)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
fi_p6 <- update(fi_p5,correlation=corARMA(p=6,form=~time|tsp))
anova(fi_p5,fi_p6)

fi_p5q1 <- update(fi_p5,correlation=corARMA(p=5,q=1,form=~time|tsp))
anova(fi_p5,fi_p5q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(fi_p5))


########################
# Plot results
########################

# First plot the raw data points for the tsp corridors
par(mfrow=c(1,1))
plot(data$time[1:192],data$fi_l[1:192],
     ylim=c(0,120),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add in the points for the control
points(data$time[193:384],data$fi_l[193:384],
       col="lightblue",
       pch=20)

# Add x-axis year labels
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Label the deployment of tsp
abline(v=89.5,lty=2)

# Calculate the seasonality offsets
offset <- mean(data$m5) * fi_p5$coef[9]
         +mean(data$m8) * fi_p5$coef[10]

# Plot the first line segment of the treated group
segments(1,
         fi_p5$coef[1] + fi_p5$coef[3] + fi_p5$coef[2] 
         + fi_p5$coef[4] + offset,
         89,
         fi_p5$coef[1] + fi_p5$coef[3] 
         + (fi_p5$coef[2] + fi_p5$coef[4])*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment of the treated group
segments(90,
         fi_p5$coef[1] + fi_p5$coef[3] 
         + (fi_p5$coef[2] + fi_p5$coef[4])*90 
         + fi_p5$coef[5] + fi_p5$coef[6] 
         + fi_p5$coef[7] + fi_p5$coef[8] + offset,
         192,
         fi_p5$coef[1] + fi_p5$coef[3] 
         + (fi_p5$coef[2] + fi_p5$coef[4])*192
         + fi_p5$coef[5] + fi_p5$coef[7]
         + (fi_p5$coef[6] + fi_p5$coef[8])*103 + offset,
         lty=1, lwd=2, col='red')

# Add the counterfactual for the treated group
segments(90, fi_p5$coef[1] + fi_p5$coef[2]*90 
         +fi_p5$coef[3] + fi_p5$coef[4]*90
         + fi_p5$coef[5] + fi_p5$coef[6] + offset,
         192, fi_p5$coef[1] + fi_p5$coef[2]*192 
         + fi_p5$coef[3] + fi_p5$coef[4]*90 
         + fi_p5$coef[5] + fi_p5$coef[6]*103 + offset,
         lty=2,col='red',lwd=2)

# Plot the first line segment of the control group
segments(1,
         fi_p5$coef[1] + fi_p5$coef[2] + offset,
         89,
         fi_p5$coef[1] + fi_p5$coef[2]*89 + offset,
         lty=1, lwd=2, col='blue')

# Plot the second line segment of the control group
segments(90,
         fi_p5$coef[1] + fi_p5$coef[2]*90 + 
           fi_p5$coef[5] + fi_p5$coef[6] + offset,
         192,
         fi_p5$coef[1] + fi_p5$coef[2]*192 + 
           fi_p5$coef[5] + fi_p5$coef[6]*103 + offset,
         lty=1, lwd=2, col='blue')

# Add the counterfactual for the control group
segments(1, fi_p5$coef[1] + fi_p5$coef[2] + offset,
         192,fi_p5$coef[1] + fi_p5$coef[2]*192 + offset,
         lty=2,col='blue',lwd=2)

# Add in a legend
legend("topright", inset=.02, legend=c("tsp","control"), col=c("red","blue"),pch=20)


##############################################
# Predict absolute and relative changes
##############################################

# Predicted crashes after tsp implementation

for (i in 1:103){
  pred <- fitted(fi_p5)[89+i]
  print(unname(pred))
}

# Estimate the counterfactual for the same period

for (i in 1:103){
  cfac <- fi_p5$coef[1] + fi_p5$coef[2]*(89+i) +
          fi_p5$coef[3] + fi_p5$coef[4]*90 +
          fi_p5$coef[5] + fi_p5$coef[6]*i
  print(unname(cfac))
}

#END



########################
### PED CRASH CONTROLLED
########################

ped_ols <- lm(ped_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m6+m7+m8+m9, data=data)
summary(ped_ols)

# Durbin-watson test, 12 time periods
dwt(ped_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for seripedy correlated errors
plot(data$time[1:192],
     residuals(ped_ols)[1:192],
     type='o',
     pch=16,
     xlab='Time',
     ylab='OLS Residuals',
     col="red")
abline(h=0,lty=2)

# Plot ACF and PACF
# Set plotting to two records on one page
par(mfrow=c(1,2))

# Produce plots
acf(residuals(ped_ols))
acf(residuals(ped_ols),type='partial')
# Note decay in ACF, no significant spike in PACF, model p=0

########################
# Run the final model
########################

# Fit the GLS regression model
ped <- gls(ped_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m6+m7+m8+m9,
              data=data)
summary(ped)
confint(ped)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
ped_p1 <- update(ped,correlation=corARMA(p=1,form=~time|tsp))
anova(ped,ped_p1)

ped_q1 <- update(ped,correlation=corARMA(q=1,form=~time|tsp))
anova(ped,ped_q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(ped))


########################
# Plot results
#########################

# First plot the raw data points for the tsp corridors
par(mfrow=c(1,1))
plot(data$time[1:192],data$ped_l[1:192],
     ylim=c(0,20),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add in the points for the control
points(data$time[193:384],data$ped_l[193:384],
       col="lightblue",
       pch=20)

# Add x-axis year labels
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Label the deployment of tsp
abline(v=89.5,lty=2)

# Calculate the seasonality offsets
offset <- (mean(data$m6)*ped$coef[9] 
           + mean(data$m7)*ped$coef[10]
           + mean(data$m8)*ped$coef[11]
           + mean(data$m9)*ped$coef[12])

# Plot the first line segment of the treated group
segments(1,
         ped$coef[1] + ped$coef[3] + ped$coef[2] 
         + ped$coef[4] + offset,
         89,
         ped$coef[1] + ped$coef[3] 
         + (ped$coef[2] + ped$coef[4])*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment of the treated group
segments(90,
         ped$coef[1] + ped$coef[3] 
         + (ped$coef[2] + ped$coef[4])*90 
         + ped$coef[5] + ped$coef[6] 
         + ped$coef[7] + ped$coef[8] + offset,
         192,
         ped$coef[1] + ped$coef[3] 
         + (ped$coef[2] + ped$coef[4])*192
         + ped$coef[5] + ped$coef[7]
         + (ped$coef[6] + ped$coef[8])*103 + offset,
         lty=1, lwd=2, col='red')

# Add the counterfactual for the treated group
segments(90, ped$coef[1] + ped$coef[2]*90 
         +ped$coef[3] + ped$coef[4]*90
         + ped$coef[5] + ped$coef[6] + offset,
         192, ped$coef[1] + ped$coef[2]*192 
         + ped$coef[3] + ped$coef[4]*90 
         + ped$coef[5] + ped$coef[6]*103 + offset,
         lty=2,col='red',lwd=2)

# Plot the first line segment of the control group
segments(1,
         ped$coef[1] + ped$coef[2] + offset,
         89,
         ped$coef[1] + ped$coef[2]*89 + offset,
         lty=1, lwd=2, col='blue')

# Plot the second line segment of the control group
segments(90,
         ped$coef[1] + ped$coef[2]*90 + 
           ped$coef[5] + ped$coef[6] + offset,
         192,
         ped$coef[1] + ped$coef[2]*192 + 
           ped$coef[5] + ped$coef[6]*103 + offset,
         lty=1, lwd=2, col='blue')

# Add the counterfactual for the control group
segments(1, ped$coef[1] + ped$coef[2] + offset,
         192,ped$coef[1] + ped$coef[2]*192 + offset,
         lty=2,col='blue',lwd=2)

# Add in a legend
legend("topright", inset=.02, legend=c("tsp","control"), col=c("red","blue"),pch=20)


##############################################
# Predict absolute and relative changes
##############################################

# Predicted crashes after tsp implementation

for (i in 1:103){
  pred <- fitted(ped)[89+i]
  print(unname(pred))
}

# Estimate the counterfactual for the same period

for (i in 1:103){
  cfac <- ped$coef[1] + ped$coef[2]*(89+i) +
          ped$coef[3] + ped$coef[4]*90 +
          ped$coef[5] + ped$coef[6]*i
  print(unname(cfac))
}

#END


########################
### BIKE CRASH CONTROLLED
########################

bike_ols <- lm(bike_l~time+tsp+tsp_time
              +level+trend+tsp_level+tsp_trend
              +m7+m8+m9, data=data)
summary(bike_ols)

# Durbin-watson test, 12 time periods
dwt(bike_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for seribikey correlated errors
plot(data$time[1:192],
     residuals(bike_ols)[1:192],
     type='o',
     pch=16,
     xlab='Time',
     ylab='OLS Residuals',
     col="red")
abline(h=0,lty=2)

# Plot ACF and PACF
# Set plotting to two records on one page
par(mfrow=c(1,2))

# Produce plots
acf(residuals(bike_ols))
acf(residuals(bike_ols),type='partial')
# Note decay in ACF, significant spike at 8 in PACF, model p=8

########################
# Run the final model
########################

# Fit the GLS regression model
bike_p8 <- gls(bike_l~time+tsp+tsp_time
           +level+trend+tsp_level+tsp_trend
           +m7+m8+m9,
           data=data,
           correlation=corARMA(p=8,form=~time|tsp),
           method="ML")
summary(bike_p8)
confint(bike_p8)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
bike_p9 <- update(bike_p8,correlation=corARMA(p=9,form=~time|tsp))
anova(bike_p8,bike_p9)

bike_p8q1 <- update(bike_p8,correlation=corARMA(p=8,q=1,form=~time|tsp))
anova(bike_p8,bike_p8q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(bike_p8))


#########################
# Plot results
#########################

# First plot the raw data points for the tsp corridors
par(mfrow=c(1,1))
plot(data$time[1:192],data$bike_l[1:192],
     ylim=c(0,10),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add in the points for the control
points(data$time[193:384],data$bike_l[193:384],
       col="lightblue",
       pch=20)

# Add x-axis year labels
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Label the deployment of tsp
abline(v=89.5,lty=2)

# Calculate the seasonality offsets
offset <- (mean(data$m7)*bike_p8$coef[9] 
           + mean(data$m8)*bike_p8$coef[10]
           + mean(data$m9)*bike_p8$coef[11])

# Plot the first line segment of the treated group
segments(1,
         bike_p8$coef[1] + bike_p8$coef[3] + bike_p8$coef[2] 
         + bike_p8$coef[4] + offset,
         89,
         bike_p8$coef[1] + bike_p8$coef[3] 
         + (bike_p8$coef[2] + bike_p8$coef[4])*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment of the treated group
segments(90,
         bike_p8$coef[1] + bike_p8$coef[3] 
         + (bike_p8$coef[2] + bike_p8$coef[4])*90 
         + bike_p8$coef[5] + bike_p8$coef[6] 
         + bike_p8$coef[7] + bike_p8$coef[8] + offset,
         192,
         bike_p8$coef[1] + bike_p8$coef[3] 
         + (bike_p8$coef[2] + bike_p8$coef[4])*192
         + bike_p8$coef[5] + bike_p8$coef[7]
         + (bike_p8$coef[6] + bike_p8$coef[8])*103 + offset,
         lty=1, lwd=2, col='red')

# Add the counterfactual for the treated group
segments(90, bike_p8$coef[1] + bike_p8$coef[2]*90 
         +bike_p8$coef[3] + bike_p8$coef[4]*90
         + bike_p8$coef[5] + bike_p8$coef[6] + offset,
         192, bike_p8$coef[1] + bike_p8$coef[2]*192 
         + bike_p8$coef[3] + bike_p8$coef[4]*90 
         + bike_p8$coef[5] + bike_p8$coef[6]*103 + offset,
         lty=2,col='red',lwd=2)

# Plot the first line segment of the control group
segments(1,
         bike_p8$coef[1] + bike_p8$coef[2] + offset,
         89,
         bike_p8$coef[1] + bike_p8$coef[2]*89 + offset,
         lty=1, lwd=2, col='blue')

# Plot the second line segment of the control group
segments(90,
         bike_p8$coef[1] + bike_p8$coef[2]*90 + 
           bike_p8$coef[5] + bike_p8$coef[6] + offset,
         192,
         bike_p8$coef[1] + bike_p8$coef[2]*192 + 
           bike_p8$coef[5] + bike_p8$coef[6]*103 + offset,
         lty=1, lwd=2, col='blue')

# Add the counterfactual for the control group
segments(1, bike_p8$coef[1] + bike_p8$coef[2] + offset,
         192,bike_p8$coef[1] + bike_p8$coef[2]*192 + offset,
         lty=2,col='blue',lwd=2)

# Add in a legend
legend("topright", inset=.02, legend=c("tsp","control"), col=c("red","blue"),pch=20)


##############################################
# Predict absolute and relative changes
##############################################

# Predicted crashes after tsp implementation

for (i in 1:103){
  pred <- fitted(bike_p8)[89+i]
  print(unname(pred))
}

# Estimate the counterfactual for the same period

for (i in 1:103){
  cfac <- bike_p8$coef[1] + bike_p8$coef[2]*(89+i) +
          bike_p8$coef[3] + bike_p8$coef[4]*90 +
          bike_p8$coef[5] + bike_p8$coef[6]*i
  print(unname(cfac))
}

#END
