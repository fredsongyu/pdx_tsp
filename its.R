library(nlme)
library(car)
library(caret)
library(tidyverse)

data <- read.delim("./data_102018.txt")
data_t <- data[1:192,]

my_data <- data_t[, c(5,6,7,8,9,16,17,18,19)]
cor(my_data,method = c("pearson"))

########################
### ALL CRASH SINGLE
########################

all_ols <- lm(tot_l~time+level+trend
                +m9, data=data_t)
summary(all_ols)

# Durbin-watson test, 12 time periods
dwt(all_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serially correlated errors
plot(data_t$time,
     residuals(all_ols),
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

########################
# Run the final model
########################

# Fit the GLS regression model
all_p1 <- gls(tot_l~time+level+trend+m9,
                 data=data_t,
                 correlation=corARMA(p=1,form=~time),
                 method="ML")

summary(all_p1)
confint(all_p1)

vif(all_p1)


########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
all_p2 <- update(all_p1,correlation=corARMA(p=2,form=~time))
anova(all_p1,all_p2)

all_p1q1 <- update(all_p1,correlation=corARMA(p=1,q=1,form=~time))
anova(all_p1,all_p1q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(all_p1))

########################
# Plot results 
########################

# Produce the plot, first plotting the raw data points
par(mfrow=c(1,1))
plot(data_t$time,data_t$tot_l,
     ylim=c(0,250),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add x axis with dates
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Add line indicating weather pattern change
abline(v=89.5,lty="dotted")

# Calculate the offset for September
offset <- mean(data_t$m9) * all_p1$coef[5]

# Plot the first line segment
segments(1,
         all_p1$coef[1] + all_p1$coef[2] + offset,
         89,
         all_p1$coef[1] + all_p1$coef[2]*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment
segments(90,
         all_p1$coef[1] + all_p1$coef[2]*90 + 
           all_p1$coef[3] + all_p1$coef[4] + offset,
         192,
         all_p1$coef[1] + all_p1$coef[2]*192 + 
           all_p1$coef[3] + all_p1$coef[4]*103 + offset,
         lty=1, lwd=2, col='red')

# Plot the counterfactual
segments(1,
         all_p1$coef[1]+all_p1$coef[2] + offset,
         192,
         all_p1$coef[1]+all_p1$coef[2]*192 + offset,
         lty=2, lwd=2, col='red')

# END



########################
### PDO CRASH SINGLE
########################

pdo_ols <- lm(pdo_l~time+level+trend
              +m9, data=data_t)
summary(pdo_ols)

# Durbin-watson test, 12 time periods
dwt(pdo_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serially correlated errors
plot(data_t$time,
     residuals(pdo_ols),
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

########################
# Run the final model
########################

# Fit the GLS regression model
pdo_p3 <- gls(pdo_l~time+level+trend+m9,
              data=data_t,
              correlation=corARMA(p=3,form=~time),
              method="ML")
summary(pdo_p3)
confint(pdo_p3)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
pdo_p4 <- update(pdo_p3,correlation=corARMA(p=4,form=~time))
anova(pdo_p3,pdo_p4)

pdo_p3q1 <- update(pdo_p3,correlation=corARMA(p=3,q=1,form=~time))
anova(pdo_p3,pdo_p3q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(pdo_p3))

########################
# Plot results 
########################

# Produce the plot, first plotting the raw data points
par(mfrow=c(1,1))
plot(data_t$time,data_t$pdo_l,
     ylim=c(0,150),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add x axis with dates
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Add line indicating weather pattern change
abline(v=89.5,lty="dotted")

# Calculate the offset for September
offset <- mean(data_t$m9) * pdo_p3$coef[5]

# Plot the first line segment
segments(1,
         pdo_p3$coef[1] + pdo_p3$coef[2] + offset,
         89,
         pdo_p3$coef[1] + pdo_p3$coef[2]*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment
segments(90,
         pdo_p3$coef[1] + pdo_p3$coef[2]*90 + 
           pdo_p3$coef[3] + pdo_p3$coef[4] + offset,
         192,
         pdo_p3$coef[1] + pdo_p3$coef[2]*192 + 
           pdo_p3$coef[3] + pdo_p3$coef[4]*103 + offset,
         lty=1, lwd=2, col='red')

# Plot the counterfactual
segments(1,
         pdo_p3$coef[1]+pdo_p3$coef[2] + offset,
         192,
         pdo_p3$coef[1]+pdo_p3$coef[2]*192 + offset,
         lty=2, lwd=2, col='red')

# END



########################
### FI CRASH SINGLE
########################

fi_ols <- lm(fi_l~time+level+trend
              +m2, data=data_t)
summary(fi_ols)

# Durbin-watson test, 12 time periods
dwt(fi_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serially correlated errors
plot(data_t$time,
     residuals(fi_ols),
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

########################
# Run the final model
########################

# Fit the GLS regression model
fi_p1 <- gls(fi_l~time+level+trend+m2,
              data=data_t,
              correlation=corARMA(p=1,form=~time),
              method="ML")
summary(fi_p1)
confint(fi_p1)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
fi_p2 <- update(fi_p1,correlation=corARMA(p=2,form=~time))
anova(fi_p1,fi_p2)

fi_p1q1 <- update(fi_p1,correlation=corARMA(p=1,q=1,form=~time))
anova(fi_p1,fi_p1q1)

# Residual plot
# Null Hypo: the residuals of a correctly specified model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(fi_p1))

########################
# Plot results 
########################

# Produce the plot, first plotting the raw data points
par(mfrow=c(1,1))
plot(data_t$time,data_t$fi_l,
     ylim=c(0,120),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add x axis with dates
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Add line indicating weather pattern change
abline(v=89.5,lty="dotted")

# Calculate the offset for September
offset <- mean(data_t$m2) * fi_p1$coef[5]

# Plot the first line segment
segments(1,
         fi_p1$coef[1] + fi_p1$coef[2] + offset,
         89,
         fi_p1$coef[1] + fi_p1$coef[2]*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment
segments(90,
         fi_p1$coef[1] + fi_p1$coef[2]*90 + 
           fi_p1$coef[3] + fi_p1$coef[4] + offset,
         192,
         fi_p1$coef[1] + fi_p1$coef[2]*192 + 
           fi_p1$coef[3] + fi_p1$coef[4]*103 + offset,
         lty=1, lwd=2, col='red')

# Plot the counterfactual
segments(1,
         fi_p1$coef[1]+fi_p1$coef[2] + offset,
         192,
         fi_p1$coef[1]+fi_p1$coef[2]*192 + offset,
         lty=2, lwd=2, col='red')

# END



########################
### PED CRASH SINGLE
########################

ped_ols <- lm(ped_l~time+level+trend
             +m4+m6+m7+m8, data=data_t)
summary(ped_ols)

# Durbin-watson test, 12 time periods
dwt(ped_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serially correlated errors
plot(data_t$time,
     residuals(ped_ols),
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

########################
# Run the final model
########################

# pedt the GLS regression model
ped <- gls(ped_l~time+level+trend
             +m4+m6+m7+m8,
             data=data_t)
summary(ped)
confint(ped)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
ped_p1 <- update(ped,correlation=corARMA(p=1,form=~time))
anova(ped,ped_p1)

ped_q1 <- update(ped,correlation=corARMA(q=1,form=~time))
anova(ped,ped_q1)

# Residual plot
# Null Hypo: the residuals of a correctly specipeded model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(ped))

########################
# Plot results 
########################

# Produce the plot, pedrst plotting the raw data points
par(mfrow=c(1,1))
plot(data_t$time,data_t$ped_l,
     ylim=c(0,20),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add x axis with dates
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Add line indicating weather pattern change
abline(v=89.5,lty="dotted")

# Calculate the offset for September
offset <- (mean(data_t$m4)*ped$coef[5]
           +mean(data_t$m6)*ped$coef[6]
           +mean(data_t$m7)*ped$coef[7]
           +mean(data_t$m8)*ped$coef[8])

# Plot the pedrst line segment
segments(1,
         ped$coef[1] + ped$coef[2] + offset,
         89,
         ped$coef[1] + ped$coef[2]*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment
segments(90,
         ped$coef[1] + ped$coef[2]*90 + 
           ped$coef[3] + ped$coef[4] + offset,
         192,
         ped$coef[1] + ped$coef[2]*192 + 
           ped$coef[3] + ped$coef[4]*103 + offset,
         lty=1, lwd=2, col='red')

# Plot the counterfactual
segments(1,
         ped$coef[1]+ped$coef[2] + offset,
         192,
         ped$coef[1]+ped$coef[2]*192 + offset,
         lty=2, lwd=2, col='red')

# END



########################
### BIKE CRASH SINGLE
########################

bike_ols <- lm(bike_l~time+level+trend
              +m5+m7+m8+m9+m10, data=data_t)
summary(bike_ols)

# Durbin-watson test, 12 time periods
dwt(bike_ols,max.lag=12,alternative="two.sided")

# Graph the residuals from the OLS regression to check for serially correlated errors
plot(data_t$time,
     residuals(bike_ols),
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

########################
# Run the final model
########################

# biket the GLS regression model
bike <- gls(bike_l~time+level+trend
            +m5+m7+m8+m9+m10,
           data=data_t)
summary(bike)
confint(bike)

########################
# Diagnostic tests
########################

# Likelihood-ratio tests to check AR process
bike_p1 <- update(bike,correlation=corARMA(p=1,form=~time))
anova(bike,bike_p1)

bike_q1 <- update(bike,correlation=corARMA(q=1,form=~time))
anova(bike,bike_q1)

# Residual plot
# Null Hypo: the residuals of a correctly specibikeed model are independently distributed--the residuals are white noise
par(mfrow=c(1,1))
qqPlot(residuals(bike))

########################
# Plot results 
########################

# Produce the plot, bikerst plotting the raw data points
par(mfrow=c(1,1))
plot(data_t$time,data_t$bike_l,
     ylim=c(0,10),
     xlim=c(0,200),
     ylab="crashes per 100 lane-miles",
     xlab="month",
     pch=20,
     col="pink",
     xaxt="n")

# Add x axis with dates
axis(1)
axis(3,at="89.5",label="tsp intervention")

# Add line indicating weather pattern change
abline(v=89.5,lty="dotted")

# Calculate the offset for September
offset <- (mean(data_t$m5)*bike$coef[5]
           +mean(data_t$m7)*bike$coef[6]
           +mean(data_t$m8)*bike$coef[7]
           +mean(data_t$m9)*bike$coef[8]
           +mean(data_t$m10)*bike$coef[9])

# Plot the bikerst line segment
segments(1,
         bike$coef[1] + bike$coef[2] + offset,
         89,
         bike$coef[1] + bike$coef[2]*89 + offset,
         lty=1, lwd=2, col='red')

# Plot the second line segment
segments(90,
         bike$coef[1] + bike$coef[2]*90 + 
           bike$coef[3] + bike$coef[4] + offset,
         192,
         bike$coef[1] + bike$coef[2]*192 + 
           bike$coef[3] + bike$coef[4]*103 + offset,
         lty=1, lwd=2, col='red')

# Plot the counterfactual
segments(1,
         bike$coef[1]+bike$coef[2] + offset,
         192,
         bike$coef[1]+bike$coef[2]*192 + offset,
         lty=2, lwd=2, col='red')

# END