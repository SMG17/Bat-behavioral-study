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
gam1 <- gam(Grooming_TotalDur_after ~ Group 
                + s(Day, by = Group)
                + s(Bat_ID, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
                data = Data_Grooming_TotalDur_after)
summary(gam1)
# Inspect model residuals:
check_resid(gam1)
gam.check(gam1)
plot(gam1, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Bat_ID random intercept:
gam2 <- gam(Grooming_TotalDur_after ~ Group 
                 + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
                 data = Data_Grooming_TotalDur_after)
summary(gam2)
# Inspect model residuals:
check_resid(gam2)
gam.check(gam2)
plot(gam2, all.terms=F, residuals=T, pch=20) 

# Run model with Day by Group smoothing term and Bat_ID random intercept:
gam3 <- gam(Grooming_TotalDur_after ~ s(Day, by = Group) 
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after)
summary(gam3)
# Inspect model residuals:
check_resid(gam3)
gam.check(gam3)
plot(gam3, all.terms=F, residuals=T, pch=20)

# Run model with Group fixed effect and Day by Group smoothing term:
gam4 <- gam(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group), 
               data = Data_Grooming_TotalDur_after)
summary(gam4)
# Inspect model residuals:
check_resid(gam4)
gam.check(gam4)
plot(gam4, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gam5 <- gam(Grooming_TotalDur_after ~ OFGroup 
               + s(Day) + s(Day, by = OFGroup) # interaction term for ordered factor
               + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
               data = Data_Grooming_TotalDur_after)
summary(gam5)
# Inspect model residuals:
check_resid(gam5)
gam.check(gam5)
plot(gam5, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gam1, gam2, gam3, gam4, gam5)
anova.gam(gam1, gam2, gam3, gam4, gam5, test="F")

# (5.1.5) run GAMM models with gamma distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gam11 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group)
             + s(Bat_ID, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gam11)
# Inspect model residuals:
check_resid(gam11)
gam.check(gam11)
plot(gam11, all.terms=F, residuals=T, pch=20) 

# Run model with Group fixed effect and Bat_ID random intercept:
gam12 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gam12)
# Inspect model residuals:
check_resid(gam12)
gam.check(gam12)
plot(gam12, all.terms=F, residuals=T, pch=20) 

# Run model with Day by Group smoothing term and Bat_ID random intercept:
gam13 <- gam(Grooming_TotalDur_after ~ s(Day, by = Group) 
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gam13)
# Inspect model residuals:
check_resid(gam13)
gam.check(gam13)
plot(gam13, all.terms=F, residuals=T, pch=20)

# Run model with Group fixed effect and Day by Group smoothing term:
gam14 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group), # interaction term for factor (Group)
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gam14)
# Inspect model residuals:
check_resid(gam14)
gam.check(gam14)
plot(gam14, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gam15 <- gam(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup) # interaction term for ordered factor
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = Gamma())
summary(gam15)
# Inspect model residuals:
check_resid(gam15)
gam.check(gam15)
plot(gam15, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gam11, gam12, gam13, gam14, gam15)
anova.gam(gam11, gam12, gam13, gam14, gam15, test="F")

# (5.1.6) run GAMM models with poisson distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gam21 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group)
             + s(Bat_ID, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gam21)
anova(gam21)

gamm21 <- gamm(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group),
              random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
              data = Data_Grooming_TotalDur_after,
              family = poisson())
summary(gamm21$gam)
anova(gamm21)

# Inspect distribution of model residuals:
check_resid(gam21)
gam.check(gam21)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gam21)^2)                                   # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-length(coef(gam21))     # estimated resid df (N-p)
resid.ssq/resid.df    # 435.7526                                       # ratio should be approx 1

plot(gam21, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Bat_ID random intercept:
gam22 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Bat_ID, Day, bs = "re"), # # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gam22)
# Inspect model residuals:
check_resid(gam22)
gam.check(gam22)
plot(gam22, all.terms=F, residuals=T, pch=20) 

# Run model with Day by Group smoothing term and Bat_ID random intercept:
gam23 <- gam(Grooming_TotalDur_after ~ s(Day, by = Group) 
             + s(Bat_ID, Day, bs = "re"), # # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gam23)
# Inspect model residuals:
check_resid(gam23)
gam.check(gam23)
plot(gam23, all.terms=F, residuals=T, pch=20)

# Run model with Group fixed effect and Day by Group smoothing term:
gam24 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group), # interaction term for factor (Group)
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gam24)
# Inspect model residuals:
check_resid(gam24)
gam.check(gam24)
plot(gam24, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gam25 <- gam(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup) # interaction term for ordered factor
             + s(Bat_ID, Day, bs = "re"), 
             data = Data_Grooming_TotalDur_after,
             family = poisson())
summary(gam25)
# Inspect model residuals:
check_resid(gam25)
gam.check(gam25)
plot(gam25, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gam21, gam22, gam23, gam24, gam25)
anova.gam(gam21, gam22, gam23, gam24, gam25, test = "Chisq")

# (5.1.7) run GAMM models with quasipoisson distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gam31 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group)
             + s(Bat_ID, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = c("quasipoisson"))
summary(gam31)
anova(gam31)

gamm31 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group),
               random = list(Bat_ID =~ 1), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = c("quasipoisson"))
summary(gam31)
anova(gam31)

# Inspect distribution of model residuals:
check_resid(gam31)
gam.check(gam31)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gam31, type = "pearson")^2)                 # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-length(coef(gam31))     # estimated resid df (N-p)
resid.ssq/resid.df  # 446.8563                                         # ratio should be approx 1

plot(gam31, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Bat_ID random intercept:
gam32 <- gam(Grooming_TotalDur_after ~ Group 
              + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
              data = Data_Grooming_TotalDur_after,
              family = c("quasipoisson"))
summary(gam32)
# Inspect model residuals:
check_resid(gam32)
gam.check(gam32)
plot(gam32, all.terms=F, residuals=T, pch=20) 

# Run model with Day by Group smoothing term and Bat_ID random intercept:
gam33 <- gam(Grooming_TotalDur_after ~ s(Day, by = Group) 
              + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
              data = Data_Grooming_TotalDur_after,
              family = c("quasipoisson"))
summary(gam33)
# Inspect model residuals:
check_resid(gam33)
gam.check(gam33)
plot(gam33, all.terms=F, residuals=T, pch=20)

# Run model with Group fixed effect and Day by Group smoothing term:
gam34 <- gam(Grooming_TotalDur_after ~ Group 
              + s(Day, by = Group), # interaction term for factor (Group)
              data = Data_Grooming_TotalDur_after,
              family = c("quasipoisson"))
summary(gam34)
# Inspect model residuals:
check_resid(gam34)
gam.check(gam34)
plot(gam34, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gam35 <- gam(Grooming_TotalDur_after ~ OFGroup 
              + s(Day) + s(Day, by = OFGroup) # interaction term for ordered factor
              + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
              data = Data_Grooming_TotalDur_after,
              family = c("quasipoisson"))
summary(gam35)
# Inspect model residuals:
check_resid(gam35)
gam.check(gam35)
plot(gam35, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gam31, gam32, gam33, gam34, gam35)
anova.gam(gam31, gam32, gam33, gam34, gam35, test = "Chisq")

# (5.1.8) run GAMM models with negative binomial distribution:
# Run model with Group fixed effect, Day by Group smoothing term, and Bat_ID random intercept:
gam41 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group)
             + s(Bat_ID, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = nb())
summary(gam41)
anova(gam41)

gamm41 <- gamm(Grooming_TotalDur_after ~ Group 
               + s(Day, by = Group),
               random = list(Bat_ID =~ 1, Bat_ID =~ nf), # Here we are including two random effects, one for just the intercept (year=~1) and another for random slope and intercept for each level of nf (year=~nf).
               data = Data_Grooming_TotalDur_after,
               family = nb())
summary(gamm41)
anova(gamm41)

# Inspect distribution of model residuals:
check_resid(gam41)
gam.check(gam41)
# Inspect dispersion of model residuals:
resid.ssq <- sum(residuals(gam41, type = "pearson")^2)                     # sum of squares of resids
resid.df <- nrow(Data_Grooming_TotalDur_after)-(length(coef(gam41))+1)     # estimated resid df (N-p)
resid.ssq/resid.df     # 0.8774                                            # ratio should be approx 1

plot(gam41, all.terms=F, shade=TRUE, scale=0, residuals=T, pch=20)
abline(h=0)

# Run model with Group fixed effect and Bat_ID random intercept:
gam42 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = nb())
summary(gam42)
# Inspect model residuals:
check_resid(gam42)
gam.check(gam42)
plot(gam42, all.terms=F, residuals=T, pch=20) 

# Run model with Day by Group smoothing term and Bat_ID random intercept:
gam43 <- gam(Grooming_TotalDur_after ~ s(Day, by = Group) 
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = nb())
summary(gam43)
# Inspect model residuals:
check_resid(gam43)
gam.check(gam43)
plot(gam43, all.terms=F, residuals=T, pch=20)

# Run model with Group fixed effect and Day by Group smoothing term:
gam44 <- gam(Grooming_TotalDur_after ~ Group 
             + s(Day, by = Group), # interaction term for factor (Group)
             data = Data_Grooming_TotalDur_after,
             family = nb())
summary(gam44)
# Inspect model residuals:
check_resid(gam44)
gam.check(gam44)
plot(gam44, all.terms=F, residuals=T, pch=20)  

# Change factor to ordered factor:
Data_Grooming_TotalDur_after$OFGroup <- as.ordered(Data_Grooming_TotalDur_after$Group)
# Change contrast to treatment coding (difference curves)
contrasts(Data_Grooming_TotalDur_after$OFGroup) <- 'contr.treatment'
# Inspect contrasts:
contrasts(Data_Grooming_TotalDur_after$OFGroup)
# Run model with OFGroup ordered factor, Day smoothing term, Day by OFGroup smoothing term, and Bat_ID random intercept:
gam45 <- gam(Grooming_TotalDur_after ~ OFGroup 
             + s(Day) + s(Day, by = OFGroup) # interaction term for ordered factor
             + s(Bat_ID, Day, bs = "re"), # random intercept s(Bat_ID, bs = "re") or slope s(Bat_ID, Day, bs = "re") (could not run random smooth s(Day, Bat_ID, bs = "fs", m = 1))
             data = Data_Grooming_TotalDur_after,
             family = nb())
summary(gam45)
# Inspect model residuals:
check_resid(gam45)
gam.check(gam45)
plot(gam45, all.terms=F, residuals=T, pch=20)  

# Compare models:
AIC(gam41, gam42, gam43, gam44, gam45)
anova.gam(gam41, gam42, gam43, gam44, gam45, test = "Chisq")

# Plot model results:
plot(gam31, all.terms=TRUE, rug=FALSE)
ggplot(data = Data_Grooming_TotalDur_after, aes(x = Day, y = Grooming_TotalDur_after, colour = Group)) +  
  geom_point() +
  geom_smooth(method= "gam", formula = y~s(x), se = F)

# version 1:
plot_smooth(gam31, view="Day", plot_all="Group")

# version 2:
plot_smooth(gam31, view="Time", cond=list(Group="Adults"), rm.ranef=TRUE, rug=FALSE, col="red", ylim=c(-15,15))
plot_smooth(gam31, view="Time", cond=list(Group="Children"), rm.ranef=TRUE, rug=FALSE, col="cyan", add=TRUE)