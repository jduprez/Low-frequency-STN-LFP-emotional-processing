
# The following code has been used to perform the statistical analyses on the data presented in the paper:
# Low-frequency subthalamic neural oscillations are involved in explicit and implicit facial emotional processing - a local field potential study, (Duprez et al., XXXX)

# This scripts needs the {lme4},{car}, {multcomp} and {MuMin} packages

library('lme4')
library('MuMIn')
library('car')
library('multcomp')

## Mixed models on behavior only - all the data
        
data = read.table(paste("/Users/joanduprez/Desktop/W/Research/UR1-EA4712/LFP/Data_LFP emotions/behavior.txt"), header=TRUE, sep = "\t")

data$emotion = as.factor(data$emotion)
data$task = as.factor(data$task)

moddacc = glmer(accuracy~emotion*task+(1|n), family=binomial(link=logit), data=data, glmerControl(optimizer = "bobyqa"))
moddacc2 = glmer(accuracy~emotion*task+(emotion|n), family=binomial(link=logit), data=data, glmerControl(optimizer = "bobyqa")) # this one
moddacc3 = glmer(accuracy~emotion*task+(task|n), family=binomial(link=logit), data=data, glmerControl(optimizer = "bobyqa"))
moddacc4 = glmer(accuracy~emotion*task+(emotion|n)+(task|n), family=binomial(link=logit), data=data, glmerControl(optimizer = "bobyqa"))

Anova(moddacc2)

qqnorm(residus)
qqline(residus)

plot(fitted(moddacc2), residuals(moddacc2),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(moddacc2), residuals(moddacc2)))

r.squaredGLMM(moddacc2)

with(data, aggregate(accuracy, list(data$emotion, data$task), mean))
with(data, aggregate(accuracy, list(data$emotion, data$task), sd))


## Power analyses

dataLFP = read.table(paste("/Users/joanduprez/Desktop/W/Research/UR1-EA4712/LFP/Data_LFP emotions/statfile_trial.txt"), header=TRUE, sep = "\t")

dataLFP$emo = as.factor(dataLFP$emo)
dataLFP$task = as.factor(dataLFP$task)
dataLFP$nucleus = as.factor(dataLFP$nucleus)

## Analyses for delta


require(lme4)
# modd1 = lmer(pwr~emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd2 = lmer(pwr~emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd3 = lmer(pwr~emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd4 = lmer(pwr~emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),]) 
modd5 = lmer(pwr~emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "delta"),]) # THIS ONE !
# modd6 = lmer(pwr~emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd7 = lmer(pwr~emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd8 = lmer(pwr~emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd9 = lmer(pwr~emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "delta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

# Get mean and sd
with(dataLFP[which(dataLFP$freq == "delta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "delta"),]$emo, dataLFP[which(dataLFP$freq == "delta"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "delta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "delta"),]$emo, dataLFP[which(dataLFP$freq == "delta"),]$task), sd))

Anova(modd5)

# Graphical check of the model's assumptions
residus<-residuals(modd5)
qqnorm(residus)
qqline(residus)

plot(fitted(modd5), residuals(modd5),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd5), residuals(modd5)))

# R squared 
r.squaredGLMM(modd5)

dataLFP$inter=interaction(dataLFP$emo,dataLFP$task)
posthoc = lmer(pwr~inter +(task|n), data=dataLFP[which(dataLFP$freq == "delta"),]) # THIS ONE !
postdelta = glht(posthoc,linfct=mcp(inter="Tukey"))
summary(postdelta)


## Analyses for alpha

# modd1 = lmer(pwr~emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
# modd2 = lmer(pwr~emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
# modd3 = lmer(pwr~emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
# modd4 = lmer(pwr~emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd5 = lmer(pwr~emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),]) # THIS ONE !
# modd6 = lmer(pwr~emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
# modd7 = lmer(pwr~emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),]) 
# modd8 = lmer(pwr~emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),]) 
# modd9 = lmer(pwr~emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "alpha"),]) 

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

Anova(modd5)
residus<-residuals(modd5)
qqnorm(residus)
qqline(residus)

plot(fitted(modd5), residuals(modd5),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd5), residuals(modd5)))

r.squaredGLMM(modd5)

with(dataLFP[which(dataLFP$freq == "alpha"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "alpha"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "alpha"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "alpha"),]$task), sd))

# Posthoc analyses
dataLFP$inter = interaction(dataLFP$emo, dataLFP$task) # Create an interaction variable
posthoc = lmer(pwr~inter +(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),])

postalpha = glht(posthoc,linfct=mcp(inter="Tukey"))
summary(postalpha)

## Analyses for beta


require(lme4)
# modd1 = lmer(pwr~emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "beta"),])
# modd2 = lmer(pwr~emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "beta"),])
# modd3 = lmer(pwr~emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "beta"),])
# modd4 = lmer(pwr~emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd5 = lmer(pwr~emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "beta"),]) # THIS ONE !
# modd6 = lmer(pwr~emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "beta"),])
# modd7 = lmer(pwr~emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "beta"),]) 
# modd8 = lmer(pwr~emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "beta"),]) 
# modd9 = lmer(pwr~emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "beta"),]) 

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

with(dataLFP[which(dataLFP$freq == "beta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "beta"),]$emo, dataLFP[which(dataLFP$freq == "beta"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "beta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "beta"),]$emo, dataLFP[which(dataLFP$freq == "beta"),]$task), sd))

Anova(modd5)
residus<-residuals(modd5)
qqnorm(residus)
qqline(residus)

plot(fitted(modd5), residuals(modd5),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd5), residuals(modd5)))


r.squaredGLMM(modd5)



## Same for ITPC

dataLFP = read.table(paste("/Users/joanduprez/Desktop/W/Research/UR1-EA4712/LFP/Data_LFP emotions/statfile_itpc.txt"), header=TRUE, sep = "\t")

dataLFP$emo = as.factor(dataLFP$emo)
dataLFP$task = as.factor(dataLFP$task)
dataLFP$nucleus = as.factor(dataLFP$nucleus)

## Analyses for delta

modd1 = lmer(pwr~emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])# THIS ONE
# modd2 = lmer(pwr~emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd3 = lmer(pwr~emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd4 = lmer(pwr~emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd5 = lmer(pwr~emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd6 = lmer(pwr~emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd7 = lmer(pwr~emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd8 = lmer(pwr~emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
# modd9 = lmer(pwr~emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "delta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

with(dataLFP[which(dataLFP$freq == "delta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "delta"),]$emo, dataLFP[which(dataLFP$freq == "delta"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "delta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "delta"),]$emo, dataLFP[which(dataLFP$freq == "delta"),]$task), sd))

Anova(modd1, test.statistic = 'F')
residus<-residuals(modd1)
qqnorm(residus)
qqline(residus)

plot(fitted(modd1), residuals(modd1),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd1), residuals(modd1)))


r.squaredGLMM(modd1)

## Analyses for theta

modd1 = lmer(pwr~emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "theta"),])# THIS ONE
# modd2 = lmer(pwr~emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd3 = lmer(pwr~emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd4 = lmer(pwr~emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd5 = lmer(pwr~emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd6 = lmer(pwr~emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd7 = lmer(pwr~emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd8 = lmer(pwr~emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
# modd9 = lmer(pwr~emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "theta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

with(dataLFP[which(dataLFP$freq == "theta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "theta"),]$emo, dataLFP[which(dataLFP$freq == "theta"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "theta"),], aggregate(pwr, list(dataLFP[which(dataLFP$freq == "theta"),]$emo, dataLFP[which(dataLFP$freq == "theta"),]$task), sd))


Anova(modd1, test.statistic = 'F')
residus<-residuals(modd1)
qqnorm(residus)
qqline(residus)

plot(fitted(modd1), residuals(modd1),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd1), residuals(modd1)))


r.squaredGLMM(modd1)

## Associations between power/ITPC and behavior (accuracy)

dataLFP = read.table(paste("yourpath/statfile_itpc.txt"), header=TRUE, sep = "\t")

## Accuracy delta 1

modd1 = lmer(accuracy~pwr+emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd2 = lmer(accuracy~pwr+emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd3 = lmer(accuracy~pwr+emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),])# THIS ONE !
modd4 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),]) 
modd5 = lmer(accuracy~pwr+emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd6 = lmer(accuracy~pwr+emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd7 = lmer(accuracy~pwr+emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd8 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd9 = lmer(accuracy~pwr+emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "delta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

with(dataLFP[which(dataLFP$freq == "delta"),], aggregate(accuracy, list(dataLFP[which(dataLFP$freq == "delta"),]$emo, dataLFP[which(dataLFP$freq == "delta"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "delta"),], aggregate(accuracy, list(dataLFP[which(dataLFP$freq == "delta"),]$emo, dataLFP[which(dataLFP$freq == "delta"),]$task), sd))

Anova(modd3, test.statistic = "F")

#can use redres package to inspect assumptions as well, only works with lme4 models
#require(redres)
#plot_redres(modd3)
#plot_resqq(modd3)
residus<-residuals(modd3)
qqnorm(residus)
qqline(residus)

plot(fitted(modd3), residuals(modd3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd3), residuals(modd3)))


r.squaredGLMM(modd3)


## Accuracy delta 2

modd1 = lmer(accuracy~pwr+emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd2 = lmer(accuracy~pwr+emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd3 = lmer(accuracy~pwr+emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "alpha"),])# THIS ONE !
modd4 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "alpha"),]) 
modd5 = lmer(accuracy~pwr+emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd6 = lmer(accuracy~pwr+emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd7 = lmer(accuracy~pwr+emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd8 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "alpha"),])
modd9 = lmer(accuracy~pwr+emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "alpha"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

with(dataLFP[which(dataLFP$freq == "alpha"),], aggregate(accuracy, list(dataLFP[which(dataLFP$freq == "alpha"),]$emo, dataLFP[which(dataLFP$freq == "alpha"),]$task), mean))
with(dataLFP[which(dataLFP$freq == "alpha"),], aggregate(accuracy, list(dataLFP[which(dataLFP$freq == "alpha"),]$emo, dataLFP[which(dataLFP$freq == "alpha"),]$task), sd))

Anova(modd3, test.statistic = "F")

qqnorm(residus)
qqline(residus)

plot(fitted(modd3), residuals(modd3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd3), residuals(modd3)))

r.squaredGLMM(modd3)

## Accuracy Beta

modd1 = lmer(accuracy~pwr+emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd2 = lmer(accuracy~pwr+emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd3 = lmer(accuracy~pwr+emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "beta"),])# THIS ONE !
modd4 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd5 = lmer(accuracy~pwr+emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd6 = lmer(accuracy~pwr+emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd7 = lmer(accuracy~pwr+emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd8 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "beta"),])
modd9 = lmer(accuracy~pwr+emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "beta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC

Anova(modd3, test.statistic = "F")

qqnorm(residus)
qqline(residus)

plot(fitted(modd3), residuals(modd3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd3), residuals(modd3)))




## Same for ITPC

dataLFP = read.table(paste("yourpath/statfile_itpc.txt"), header=TRUE, sep = "\t")

dataLFP$emo = as.factor(dataLFP$emo)
dataLFP$task = as.factor(dataLFP$task)
dataLFP$nucleus = as.factor(dataLFP$nucleus)

# Delta

modd1 = lmer(accuracy~pwr+emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd2 = lmer(accuracy~pwr+emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd3 = lmer(accuracy~pwr+emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),])# THIS ONE !
modd4 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "delta"),]) 
modd5 = lmer(accuracy~pwr+emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd6 = lmer(accuracy~pwr+emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd7 = lmer(accuracy~pwr+emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd8 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "delta"),])
modd9 = lmer(accuracy~pwr+emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "delta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC


Anova(modd3, test.statistic = "F")

qqnorm(residus)
qqline(residus)

plot(fitted(modd3), residuals(modd3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd3), residuals(modd3)))

r.squaredGLMM(modd3)


# Theta

modd1 = lmer(accuracy~pwr+emo*task +(1|n), data=dataLFP[which(dataLFP$freq == "theta"),])
modd2 = lmer(accuracy~pwr+emo*task+nucleus+(1|n), data=dataLFP[which(dataLFP$freq == "theta"),])
modd3 = lmer(accuracy~pwr+emo*task +(emo|n), data=dataLFP[which(dataLFP$freq == "theta"),])# THIS ONE !
modd4 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n), data=dataLFP[which(dataLFP$freq == "theta"),]) 
modd5 = lmer(accuracy~pwr+emo*task +(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
modd6 = lmer(accuracy~pwr+emo*task+nucleus+(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
modd7 = lmer(accuracy~pwr+emo*task +(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
modd8 = lmer(accuracy~pwr+emo*task+nucleus+(emo|n)+(task|n), data=dataLFP[which(dataLFP$freq == "theta"),])
modd9 = lmer(accuracy~pwr+emo*task + nucleus +(emo|n)+(task|n)+(n|nucleus), data=dataLFP[which(dataLFP$freq == "theta"),])

anova(modd1, modd2, modd3, modd4, modd5, modd6, modd7, modd8, modd9) # Compare models and keep the one that converges with the lower AIC


Anova(modd3, test.statistic = "F")

qqnorm(residus)
qqline(residus)

plot(fitted(modd3), residuals(modd3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modd3), residuals(modd3)))

r.squaredGLMM(modd3)



