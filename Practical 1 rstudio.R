### Practical 1 ###
#install.packages("reshape")
library(ggplot2)
library(plyr)
library(rms)
library(MASS)
library(reshape)
library(lmtest)

# You might need to set your working directory here 

smart <- readRDS(file= "SMARTs_P1.rds")
class(smart) # This tells us the class/data type of the object smart (hopefully returning data.frame)
dim(smart)
sapply(smart,class) # This tells us the class of each variable in the smart data, eg, categorical='factor' or continuous='numeric'
Albumin is not a factor variable so we need to change and label

# Albumin
smart$albumin <- as.factor(smart$albumin)
levels(smart$albumin)
smart$albumin <- revalue(smart$albumin, c("1"="No", "2"="Low", "3"="High"))
levels(smart$albumin)
There are 2 variables for systolic blood pressure (SBP) because it can either done by hand or automatic. If we check, we can see there is a lot of missing data in both of the SBP variables.

table(is.na(smart$SYSTBP)==FALSE)

table(is.na(smart$SYSTH)==FALSE)

table(is.na(smart$SYSTH)==FALSE | is.na(smart$SYSTBP)==FALSE)

smart$SBP <- ifelse(is.na(smart$SYSTBP)==FALSE, smart$SYSTBP, smart$SYSTH)
table(is.na(smart$SBP))

## Part 2 Exploratory analysis

# Overview of the data
attach(smart)
dim(smart)
summary(smart)
sapply(smart,class)
describe(smart) 

# Summarise outcome
table(outcome)
round(prop.table(table(outcome))*100,0)

# Summarise categorical variables by outcome
table(SEX,outcome) # Produces a table of counts
round(prop.table(table(SEX,outcome),2)*100,0) # Produces a table of percentages
round(prop.table(table(outcome,SEX),2)*100,0) # Produces a table of percentages

table(SMOKING,outcome) 
round(prop.table(table(SMOKING,outcome),2)*100,0)
round(prop.table(table(outcome,SMOKING),2)*100,0)

# Look at distributions of continuous variables
summary(AGE)
ggplot(smart,aes(x=AGE))+geom_histogram()+facet_grid(~outcome)+theme_light()+ labs(title = "Histogram of age split by outcome")
ggplot(smart, aes(AGE, fill = outcome)) + geom_histogram()
ggplot(smart, aes(AGE, fill = as.factor(outcome))) + geom_density(alpha = 0.2)

# Question - Above we look into age, can you do the same for BMI?
#   What is the range of BMI values?
#   Plot a BMI histogram
# Plot the BMI density by outcome


## Part 3 Univariate logistic regression with categorical variables
# Lets start by producing a univariate logistic model with gender as the independant variable.

detach(smart)
# Recode SEX so that the baseline is male #
smart$SEX2 = relevel(smart$SEX, ref="Male")
attach(smart)

# logistic regression model including sex as the only independent variable
sex_mod <- glm(outcome~SEX2,family="binomial")
summary(sex_mod)
beta <- round(exp(sex_mod$coef),3) # extract the model coefficients from the model
b <- round(exp(confint.default(sex_mod)),3) # Compute the confidence intervals for the model coefficients
d <- cbind(beta,b)
d
Can you interpret the output?
  
# By modifying the code above, try fitting a logistic regression model to the SMART data with albumin as the predictor.
# Q. Interpret the output. Is albumin a significant predictor? What level of albumin is associated with the lowest risk?
  

## Part 4 Univariate logistic regression with continuous variables
# Now lets try with a continuous variable (AGE).

# Assuming age is linear
age_mod <- glm(outcome~AGE,family="binomial")
summary(age_mod)
exp(confint.default((age_mod))) # if you want confidence intervals for the model coefficients

lp40 <- predict(age_mod,data.frame(AGE=40))
risk40 <- exp(lp40)/(1+exp(lp40))*100

lp80 <- predict(age_mod,data.frame(AGE=80))
risk80 <- exp(lp80)/(1+exp(lp80))*100

matrix(c("risk40","risk80", risk40, risk80),ncol=2,byrow=F)

risk40a <- predict(age_mod,data.frame(AGE=40), type="response")*100 
# 'type="response"' gives the predicted probabilites and so multiply by 100 to get the risk as a percentage
risk40b <- 1/(1+exp(-lp40))*100 # alternative calculation
c(risk40a, risk40b)

## Part 5 Modelling a continuos variable using splines

age3_spline <- rcs(AGE,3)
age4_spline <- rcs(AGE,4)
age5_spline <- rcs(AGE,5)

lp_age1 <- predict(age_mod)

age3_mod <- glm(outcome~age3_spline,family="binomial")
age3_mod
lp_age3 <- predict(age3_mod)

age4_mod <- glm(outcome~age4_spline,family="binomial")
age4_mod
lp_age4 <- predict(age4_mod)

age5_mod <- glm(outcome~age5_spline,family="binomial")
age5_mod
lp_age5 <- predict(age5_mod)

data_part6 <- data.frame(AGE,lp_age1,lp_age3,lp_age4,lp_age5)
data_part6_m <- melt(data_part6,id.vars='AGE')
plot_part6 <- ggplot(data_part6_m,aes(AGE,value,colour=variable))+geom_line()+scale_colour_manual(labels=c("linear","3 knots","4 knots","5 knots"),values=c("gray","green","red","blue"))+theme_bw()
plot_part6 + labs(x="Age (years)",y="Linear Predictor (log odds)",color="") + theme(legend.position=c(0.2,0.8))


age_spline_check <- matrix(c(AIC(age_mod),
                             BIC(age_mod),
                             AIC(age3_mod),
                             BIC(age3_mod),
                             AIC(age4_mod),
                             BIC(age4_mod),
                             AIC(age5_mod),
                             BIC(age5_mod)), ncol=2, byrow=TRUE)

colnames(age_spline_check) <- c("AIC", "BIC")
rownames(age_spline_check) <- c("age_mod","age3_mod", "age4_mod", "age5_mod")
age_spline_check


## Part 5 Building a multivariable model

detach(smart)
smart <- subset(smart, select = c(outcome, SEX, AGE, SBP, alcohol, CHOLO, BMIO, DIABETES, CARDIAC, SMOKING, AAA))
smart <- na.omit(smart)
age3_spline <- rcs(smart$AGE,3)
attach(smart)

k10 <- qchisq(0.10,1,lower.tail=FALSE) # this give the change in AIC we consider to be significant in our stepwise selection

# Forward selection (by AIC)
empty_mod_2 <- glm(outcome~1,family="binomial")
forward_mod_2 <- stepAIC(empty_mod_2,k=k10,scope=list(upper=~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,lower=~1),direction="forward",trace=TRUE)

# Backward selection (by AIC)
full_mod_2 <- glm(outcome~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,family="binomial")
backward_mod_2 <- stepAIC(full_mod_2,k=k10,scope=list(upper=~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,lower=~1),direction="backward",trace=TRUE)

forward_mod_2
backward_mod_2

# Backward selection (by AIC) forcing SEX to be incuded in the model
backward_mod_2sex <- stepAIC(full_mod_2,k=k10,scope=list(upper=~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,lower=~SEX),direction="backward",trace=TRUE)
backward_mod_2sex
summary(backward_mod_2sex)

