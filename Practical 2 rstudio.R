### Practical 2 ###

# install.packages("ggplot2")

library(ggplot2)
library(plyr)
library(rms)
library(reshape2)
library(prodlim)
library(muhaz)
library(survival)
library(Hmisc)
library(MASS)

## Part 1

smart <- readRDS("SMARTs_P2-3.rds")
attach(smart)

class(smart)

sapply(smart,class)

dim(smart)

sapply(smart,summary)

summary(smart$TEVENT) # Summarises the time to event

table(smart$EVENT) # Tabulates number of events

table(smart$EVENT==1, round(smart$TEVENT/365,0)) # Tabulates number of events by year (TRUE=event, FALSE=censored)

## Part 2
#Survival objects and Kaplan-Meier plots

plot(KM,xlab="Time (days)",ylab="Survival")

TEmth<-Surv(TEVENT/12,EVENT)
KMmth<-survfit(TEmth~1, data = smart)

plot(KMmth,xlab="Time (months)",ylab="Survival", ylim = c(0.5,1))
 
KM <-survfit(TE~SEX, data = smart)

plot(KM,lty=1:2,xlab="Time",ylab="Survival", ylim = c(0.5,1))

#ggsurvplot(KM, data = smart, risk.table = TRUE, ylim = c(0.5,1))
 
## Part 3
#Now lets begin by fitting univariable models with categorical variables

cox_sex <- coxph(Surv(TEVENT/12,EVENT)~SEX)
summary(cox_sex)

cox_card <- coxph(Surv(TEVENT/12,EVENT)~CARDIAC)
summary(cox_card)

## Part 4
#Modelling using splines

age3_spline <- rcs(AGE,3)
age4_spline <- rcs(AGE,4)
age5_spline <- rcs(AGE,5)


# Fit a Cox model for each spline function
cox_age3 <- coxph(Surv(TEVENT/12,EVENT)~age3_spline,data=smart,ties="breslow")
cox_age3
lp_age3 <- predict(cox_age3)

cox_age4 <- coxph(Surv(TEVENT/12,EVENT)~age4_spline,data=smart,ties="breslow")
cox_age4
lp_age4 <- predict(cox_age4)

cox_age5 <- coxph(Surv(TEVENT/12,EVENT)~age5_spline,data=smart,ties="breslow")
cox_age5
lp_age5 <- predict(cox_age5)


cox_age <- coxph(Surv(TEVENT/12,EVENT)~AGE)
summary(cox_age)
lp_age1 <- predict(cox_age)


data_part6 <- data.frame(AGE,lp_age1,lp_age3,lp_age4,lp_age5)
data_part6_m <- melt(data_part6,id.vars='AGE')
plot_part6 <- ggplot(data_part6_m,aes(AGE,value,colour=variable))+geom_line()+scale_colour_manual(labels=c("linear","3 knots","4 knots","5 knots"),values=c("gray","green","red","blue"))+theme_bw()
plot_part6 + labs(x="Age (years)",y="Linear Predictor (log odds)",color="") + theme(legend.position=c(0.2,0.8))
age_spline_check <- matrix(c(AIC(cox_age),
                             BIC(cox_age),
                             AIC(cox_age3),
                             BIC(cox_age3),
                             AIC(cox_age4),
                             BIC(cox_age4),
                             AIC(cox_age5),
                             BIC(cox_age5)), ncol=2, byrow=TRUE)

colnames(age_spline_check) <- c("AIC", "BIC")
rownames(age_spline_check) <- c("age_mod","age3_mod", "age4_mod", "age5_mod")
age_spline_check

Part 5 Building a multivariable model

detach(smart)
smart <- subset(smart, select = c(outcome, TEVENT, EVENT, SEX, AGE, SBP, alcohol, CHOLO, BMIO, DIABETES, CARDIAC, SMOKING, AAA))
smart <- na.omit(smart)
age3_spline <- rcs(smart$AGE,3)
attach(smart)

k10 <- qchisq(0.10,1,lower.tail=FALSE) # this gives the change in AIC we consider to be significant in our stepwise selection

# Forward selection (by AIC)
empty_mod_2 <- coxph(Surv(TEVENT/12,EVENT)~1)
forward_mod_2 <- stepAIC(empty_mod_2,k=k10,scope=list(upper=~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,lower=~1),direction="forward",trace=TRUE)

# Backward selection (by AIC)
full_mod_2 <- coxph(Surv(TEVENT/12,EVENT)~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA)
backward_mod_2 <- stepAIC(full_mod_2,k=k10,scope=list(upper=~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,lower=~1),direction="backward",trace=TRUE)

summary(forward_mod_2)
summary(backward_mod_2)
Remember the logistic model in practical 1 choose age, AAA and diabetes as significant predictors and we forced sex to be in the final model. Here, sex was chosen as a significant predictor in both the forward and backward selection approachs.
Backward selection also included alcohol and smoking as predictors in the final model. As disscussed in practical 1, forward and backward selection do not always chose the same final model and we need to be careful when building our model.

For the stepwise selection we have been using a p-value of 0.1 (equivalent to an AIC of 2.706) so not to be too stringent. However, if we change this to 0.05 we can end with a different model. Lets try this in backward selection....

k10 <- qchisq(0.05,1,lower.tail=FALSE)

full_mod_2 <- coxph(Surv(TEVENT/12,EVENT)~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA)
backward_mod_2a <- stepAIC(full_mod_2,k=k10,scope=list(upper=~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,lower=~1),direction="backward",trace=TRUE)


