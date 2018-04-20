# Clear variables from the memory
closeAllConnections()
rm(list=ls())

############################################################
# (1) Set environment and load dataset and libraries
############################################################
## (1.1) Set environment
readPath <- setwd( "C:/Users/sarah/OneDrive - Temple University/PhD Temple 2015-2020/PhD Research/Video analysis/Behavioral data")
readPath <- setwd( "C:/Users/Marianne/OneDrive - Temple University/PhD Temple 2015-2020/PhD Research/Video analysis/Behavioral data")

## (1.2) Load .csv dataset
Data <- read.csv("Behavioral_data_by_arousal_2018_04_14.csv")
Data <- read.csv("Individual_behavioral_data_2018_04_13_FirstHalf.csv")
Data <- read.csv("Individual_behavioral_data_2018_04_13_SecondHalf.csv")
Data <- read.csv("Behavioral_data_by_arousal_without_zero_Grooming_TotalDur_after_2018_04_14.csv")

## (1.3) Load libraries
library(lme4)
library(lmerTest)
library(fitdistrplus)
library(car)
library(ggplot2)
library(mgcv)
library(itsadug)
library(AER)

############################################################
# (2) Format conversions
############################################################
## (2.1) Date and time formating
Data$First_movement <- as.POSIXct(as.character(Data$First_movement), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Arousal_start <- as.POSIXct(as.character(Data$Arousal_start), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Active_start <- as.POSIXct(as.character(Data$Active_start), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Arousal_end <- as.POSIXct(as.character(Data$Arousal_end), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Last_movement <- as.POSIXct(as.character(Data$Last_movement), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")
Data$Date_died <- as.POSIXct(as.character(Data$Date_died), format= "%Y-%m-%d %H:%M:%S", tz = "UTC")

## (2.2) Factor formating
Data$Bat_ID <- factor(Data$Bat_ID)
Data$Group <- factor(Data$Group)

## (2.3) Remove NA rows
Data <- Data[!(is.na(Data$Active_start)), ]

############################################################
# (3) Assess correlations among predictors
############################################################
## (3.1) Is ArousalDur2 correlated with Group?
# (3.1.1) Anova with Group and random effect of BatID (ArousalDur per arousal)
lmer_ArousalDur <- lmer(ArousalDur2 ~ Group + (1|Bat_ID), data = Data, na.action=na.omit)
summary(lmer_ArousalDur)
anova(lmer_ArousalDur)

# Post hoc tests
difflsmeans(lmer_ArousalDur, test.effs="Group")
boxplot(Data$ArousalDur2 ~ Data$Group, xlab = "Group", ylab = "Arousal duration (sec)")

# Normality of residuals
plot(lmer_ArousalDur)
qqnorm(resid(lmer_ArousalDur))
qqline(resid(lmer_ArousalDur))
hist(resid(lmer_ArousalDur))

# (3.1.2) Anova with Group (ArousalDur as individual means)
ArousalDur_means <- aggregate(Data$ArousalDur2, Data[,c(2,6)], FUN = mean, na.rm = TRUE) # for each arousal
colnames(ArousalDur_means)[3] <- "ArousalDur_mean"
aov_ArousalDur <- aov(ArousalDur_mean ~ Group, data = ArousalDur_means)
summary(aov_ArousalDur)
model.tables(aov_ArousalDur, "means")
boxplot(ArousalDur_means$ArousalDur_mean ~ ArousalDur_means$Group, xlab = "Group", ylab = "Mean arousal duration (sec)")

# Test assumptions
bartlett.test(ArousalDur_mean ~ Group, data = ArousalDur_means) #Bartlett test of homogeneity of variances - Significant result, therefore variances cannot be assumed to be equal
plot(aov_ArousalDur)

# Conclusion: Yes, ArousalDur is correlated with Group, so won't be included as fixed effect in following models.

############################################################
# (4) Determine probability distribution of data
############################################################
## (4.1) Probability distribution of Grooming_TotalDur_after
# (4.1.1) Remove rows with missing values
Data_Grooming_TotalDur_after <- Data[!(is.na(Data$Grooming_TotalDur_after)), ]

# (4.1.2) Assess data distribution
hist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after)
qqnorm(Data_Grooming_TotalDur_after$Grooming_TotalDur_after)
qqline(Data_Grooming_TotalDur_after$Grooming_TotalDur_after)

# (4.1.3) Fit normal probability distribution
f1 <- fitdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"norm")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"norm",para=list(mean=f1$estimate[1],sd=f1$estimate[2]))
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "norm")
f1$aic

f1 <- fitdistr(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"normal")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"norm",para=list(mean=f1$estimate[1],sd=f1$estimate[2]))
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "norm")
AIC(f1)

# (4.1.4) Fit gamma probability distribution
f2 <- fitdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"gamma")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"gamma",para=list(shape=f2$estimate[1],rate=f2$estimate[2]))
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "gamma", shape = f2$estimate[1], rate = f2$estimate[2])
f2$aic

f2 <- fitdistr(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "gamma", lower = 0.01) 
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"gamma",para=list(shape=f2$estimate[1],rate=f2$estimate[2]))
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "gamma", shape = f2$estimate[1], rate = f2$estimate[2])
AIC(f2)

# (4.1.5) Fit lognormal probability distribution
f3 <- fitdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"lnorm")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"lnorm",para=list(meanlog=f3$estimate[1],sdlog=f3$estimate[2])) #YES!
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "lnorm")
f3$aic

f3 <- fitdistr(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"lognormal")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"lnorm",para=list(meanlog=f3$estimate[1],sdlog=f3$estimate[2])) #YES!
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "lnorm")
AIC(f3)

# (4.1.6) Fit normal probability distribution on log-transformed data
hist(log(Data_Grooming_TotalDur_after$Grooming_TotalDur_after))
qqnorm(log(Data_Grooming_TotalDur_after$Grooming_TotalDur_after))
qqline(log(Data_Grooming_TotalDur_after$Grooming_TotalDur_after))

f4 <- fitdist(log(Data_Grooming_TotalDur_after$Grooming_TotalDur_after),"norm")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"norm",para=list(mean=f4$estimate[1],sd=f4$estimate[2]))
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "norm")
f4$aic

f4 <- fitdistr(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"normal")
plotdist(Data_Grooming_TotalDur_after$Grooming_TotalDur_after,"norm",para=list(mean=f4$estimate[1],sd=f4$estimate[2]))
qqp(Data_Grooming_TotalDur_after$Grooming_TotalDur_after, "norm")
AIC(f4)


############################################################
# (5) Fit generalized additive mixed modesls  
############################################################
## (5.1) GAMM with  response variable Grooming_TotalDur_after
# (5.1.1) Add 1 to grooming duration values in dataset
Data_Grooming_TotalDur_after <- Data
Data_Grooming_TotalDur_after$Grooming_TotalDur_after <- (Data_Grooming_TotalDur_after$Grooming_TotalDur_after + 1)

# (5.1.2) plot Grooming_TotalDur_after as a function of Day
plot(Grooming_TotalDur_after ~ Day, pch = 20, data = Data_Grooming_TotalDur_after)

# (5.1.3) plot Grooming_TotalDur_after as a function of Day and Group
ggplot(data = Data_Grooming_TotalDur_after, aes(x = Day, y = Grooming_TotalDur_after, colour = Group, group = Group)) +  
  geom_point(aes(col = Data_Grooming_TotalDur_after$Group)) +
  geom_smooth(se = T)

# (5.1.4) run GAMM models with normal distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm1 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group),
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after)
summary(gamm1)
# Inspect model residuals:
check_resid(gamm1)
gam.check(gamm1)
plot(gamm1, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Bat_ID random intercept:
gamm2 <- gamm(Grooming_TotalDur_after ~ Group,
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after)
summary(gamm2)
# Inspect model residuals:
check_resid(gamm2)
gam.check(gamm2)
plot(gamm2, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day by Group smoothing term:
gamm3 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group), 
               data = Data_Grooming_TotalDur_after)
summary(gamm3)
# Inspect model residuals:
check_resid(gamm3)
gam.check(gamm3)
plot(gamm3, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm4 <- gamm(Grooming_TotalDur_after ~ OFGroup 
              + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after)
summary(gamm4)
# Inspect model residuals:
check_resid(gamm4)
gam.check(gamm4)
plot(gamm4, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm1, gamm2, gamm3, gamm4)
anova.gamm(gamm1, gamm2, gamm3, gamm4, test="F")

# (5.1.5) run GAMM models with gamma distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm11 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group),
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after,
              family = Gamma())
summary(gamm11)
# Inspect model residuals:
check_resid(gamm11)
gam.check(gamm11)
plot(gamm11, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Bat_ID random intercept:
gamm12 <- gamm(Grooming_TotalDur_after ~ Group,
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = Gamma())
summary(gamm12)
# Inspect model residuals:
check_resid(gamm12)
gam.check(gamm12)
plot(gamm12, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day by Group smoothing term:
gamm13 <- gamm(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group), # interaction term for factor (Group)
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gamm13)
# Inspect model residuals:
check_resid(gamm13)
gam.check(gamm13)
plot(gamm13, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm14 <- gamm(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
             random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gamm14)
# Inspect model residuals:
check_resid(gamm14)
gam.check(gamm14)
plot(gamm14, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm11, gamm12, gamm13, gamm14)
anova.gam(gamm11, gamm12, gamm13, gamm14, test="F")

# (5.1.6) run GAMM models with poisson distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm21 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group),
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after,
              family = poisson())
summary(gamm21$gam)
anova(gamm21)

# Inspect distribution of model residuals:
check_resid(gamm21)
gam.check(gamm21)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gamm21)^2)                                   # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-length(coef(gamm21))     # estimated resid df (N-p)
resid.ssq/resid.df    # 435.7526                                       # ratio should be approx 1

plot(gamm21, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Bat_ID random intercept:
gamm22 <- gamm(Grooming_TotalDur_after ~ Group,
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = poisson())
summary(gamm22)
# Inspect model residuals:
check_resid(gamm22)
gam.check(gamm22)
plot(gamm22, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day by Group smoothing term:
gamm23 <- gamm(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group), # interaction term for factor (Group)
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gamm23)
# Inspect model residuals:
check_resid(gamm23)
gam.check(gamm23)
plot(gamm23, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm24 <- gamm(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
             random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gamm24)
# Inspect model residuals:
check_resid(gamm24)
gam.check(gamm24)
plot(gamm24, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm21, gamm22, gamm23, gamm24)
anova.gam(gamm21, gamm22, gamm23, gamm24, test = "Chisq")

# (5.1.7) run GAMM models with quasipoisson distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm31 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group),
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = c("quasipoisson"))
summary(gamm31)
anova(gamm31)

# Inspect distribution of model residuals:
check_resid(gamm31)
gam.check(gamm31)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gamm31, type = "pearson")^2)                 # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-length(coef(gamm31))     # estimated resid df (N-p)
resid.ssq/resid.df  # 446.8563                                         # ratio should be approx 1

plot(gamm31, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Bat_ID random intercept:
gamm32 <- gamm(Grooming_TotalDur_after ~ Group,
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = c("quasipoisson"))
summary(gamm32)
# Inspect model residuals:
check_resid(gamm32)
gam.check(gamm32)
plot(gamm32, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day by Group smoothing term:
gamm33 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group), # interaction term for factor (Group)
              data = Data_Grooming_TotalDur_after,
              family = c("quasipoisson"))
summary(gamm33)
# Inspect model residuals:
check_resid(gamm33)
gam.check(gamm33)
plot(gamm33, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm34 <- gamm(Grooming_TotalDur_after ~ OFGroup 
              + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after,
              family = c("quasipoisson"))
summary(gamm34)
# Inspect model residuals:
check_resid(gamm34)
gam.check(gamm34)
plot(gamm34, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm31, gamm32, gamm33, gamm34)
anova.gam(gamm31, gamm32, gamm33, gamm34, test = "Chisq")

# (5.1.8) run GAMM models with negative binomial distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gamm41 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group),
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = nb())
summary(gamm41)
anova(gamm41)

# Inspect distribution of model residuals:
check_resid(gamm41)
gam.check(gamm41)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gamm41, type = "pearson")^2)                     # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-(length(coef(gamm41))+1)     # estimated resid df (N-p)
resid.ssq/resid.df     # 0.8774                                            # ratio should be approx 1

plot(gamm41, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Bat_ID random intercept:
gamm42 <- gamm(Grooming_TotalDur_after ~ Group,
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = nb())
summary(gamm42)
# Inspect model residuals:
check_resid(gamm42)
gam.check(gamm42)
plot(gamm42, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Day by Group smoothing term:
gamm43 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group), # interaction term for factor (Group)
               data = Data_Grooming_TotalDur_after,
               family = nb())
summary(gamm43)
# Inspect model residuals:
check_resid(gamm43)
gam.check(gamm43)
plot(gamm43, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gamm44 <- gamm(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup), # interaction term for ordered factor
             random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
             data = Data_Grooming_TotalDur_after,
             family = nb())
summary(gamm44)
# Inspect model residuals:
check_resid(gamm44)
gam.check(gamm44)
plot(gamm44, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gamm41, gamm42, gamm43, gamm44)
anova.gam(gamm41, gamm42, gamm43, gamm44, test = "Chisq")

# Plot model results:
plot(gamm31, all.terms=TRUE, rug=FALSE)
ggplot(data = Data_Grooming_TotalDur_after, aes(x = Day, y = Grooming_TotalDur_after, colour = Group)) +  
  geom_point() +
  geom_smooth(method= "gam", formula = y~s(x), se = F)

# version 1:
plot_smooth(gamm31, view="Day", plot_all="Group")

# version 2:
plot_smooth(gamm31, view="Time", cond=list(Group="Adults"), rm.ranef=TRUE, rug=FALSE, col="red", ylim=c(-15,15))
plot_smooth(gamm31, view="Time", cond=list(Group="Children"), rm.ranef=TRUE, rug=FALSE, col="cyan", add=TRUE)