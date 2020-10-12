### Practical 3 - Validation ###

library(pROC)
library(DescTools)
library(Hmisc)
library(plyr)
library(boot)
library(MASS)
library(rms)

## Part 1

smart <- readRDS("SMARTs_P2-3.rds")
smart <- subset(smart, select = c(outcome, TEVENT, EVENT, SEX, AGE, SBP, alcohol, CHOLO, BMIO, DIABETES, CARDIAC, SMOKING, AAA))
smart <- na.omit(smart)
attach(smart)

class(smart)

sapply(smart,class)

## Part 2
#Fit the multivariable models selected from practical 1

age3_spline <- rcs(AGE,3)

# Fit model selected from practicals 1 and 2
mod_log <- glm(outcome~SEX+SBP+DIABETES+age3_spline,family="binomial",x=TRUE,y=TRUE)
mod_log

## Part 3 apparent discrimination

## Calculate apparent discrimination performance
# Obtain the predicted probabilities for each patient
pred_prob <- predict(mod_log,type="response")

# Obtain the c statistic / AUC
c1 <- roc(outcome~pred_prob,ci=TRUE)
c1

# Plot the ROC curve
plot(roc(outcome,pred_prob))

#Brier score
BrierScore(mod_log)

## Part 4 apparent calibration

# Obtain the linear predictor/PI for each patient
pred_LP <- predict(mod_log,type="link")

# CITL
mod_log_1 <- glm(outcome~offset(pred_LP),family="binomial")
coefficients(mod_log_1)
confint(mod_log_1)

# E/O
prop.table(table(outcome))
O <- prop.table(table(outcome))[2]
E <- mean(pred_prob)
E/O


mod_log_2 <- glm(outcome~pred_LP,family="binomial",x=TRUE,y=TRUE)
coefficients(mod_log_2)
confint(mod_log_2)
 
# Visual assessment of calibration by risk groups
# create 10 risk groups
groups <- cut(pred_prob,breaks=quantile(pred_prob, prob = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)),labels=c(1:10),include.lowest=TRUE)

# average the observed and expected probabilities of patients in each risk group 
gpdata <- cbind(smart,groups,pred_prob)
obs <- (ddply(gpdata,~groups,summarise,mean=mean(as.numeric(outcome)))[,2])
exp <- ddply(gpdata,~groups,summarise,mean=mean(pred_prob))
obsn <- table(outcome,groups)[1,] 

# CIs for scatter points
lci = pmax(0,(obs - (1.96*(((obs*(1-obs))/obsn)^.5))))
uci = pmin(1,(obs + (1.96*(((obs*(1-obs))/obsn)^.5))))

# graph a simple calibration plot over 10 risk groups
par(pty="s")
plot(obs~exp[,2],xlim=c(0,1),ylim=c(0,1),col="red",ylab="Observed",xlab="Expected")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}
h <- hist(pred_prob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(outcome))~pred_prob,span=1))
lines_data <- data.frame(pred_prob,obs_all)
lines_data2 <- lines_data[order(pred_prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.9,c("Risk groups","Reference line","95% CI","Loess"),col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")

par(pty="s")
plot(obs~exp[,2],xlim=c(0,.3),ylim=c(0,.3),col="red",ylab="Observed",xlab="Expected")
lines(c(0,1),c(0,1),lty=2)
for(i in 1:10){
  lines(c(exp[i,2],exp[i,2]),c(lci[i],uci[i]),col="green")
}
h <- hist(pred_prob, breaks=50, plot=FALSE)
for(i in 1:length(h$mids)){
  lines(c(h$mids[i],h$mids[i]),c(rep(1,length(h$mids))[i],1-((h$counts[i]/max(h$counts))/10)))
}

# Add a loess smoother to the plot
obs_all <- predict(loess((as.numeric(outcome))~pred_prob,span=1))
lines_data <- data.frame(pred_prob,obs_all)
lines_data2 <- lines_data[order(pred_prob),] 
lines(lines_data2[,1],lines_data2[,2],col="blue")
legend(0.0,0.25,c("Risk groups","Reference line","95% CI","Loess"),col=c("red","black","green","blue"),lty=c(0,2,1,1),pch=c(1,NA,NA,NA),bty="n")

## Part 5 bootstrap validation

# You can obtain c-statistic and c-slope directly from Harrell's "validate" programme within 'rms'.
mod1 <- lrm(outcome~SEX+SBP+DIABETES+age3_spline,x=TRUE,y=TRUE)
set.seed(231398)
boot_1 <- validate(mod1,method="boot",B=50)
boot_1

# Note that this gives Dxy rather than c, however Dxy = 2*(c-0.5), i.e. c=(Dxy/2)+0.5
(boot_1[1,1]+1)/2
(boot_1[1,5]+1)/2

## To get bootstrapped answers for CITL, c-statistic & c-slope it's a more manual approach:
manual_boot <- function(data,samples){
  results <- matrix(nrow = samples,ncol = 6)
  set.seed(231398)
  for (i in 1:samples) {
    samp_index <- sample(1:nrow(data), nrow(data), rep=TRUE) # create a sampling index vector
    
    bs_samp <- data[samp_index,] # index the orignal dataset using the sampling vector to give the bs sample
    model <- glm(outcome~SEX+age3_spline + DIABETES,family=binomial, data=bs_samp)  # Fit model to the bootstrap sample
    pr_bs <- predict(model,type="response") 
    lp_bs <- predict(model) # predict lp from the bootstrap model in the bs sample
    
    pr_test <- predict(model,type="response",newdata = data)  # predict probabilities from the bootstrap model in the original sample
    lp_test <- predict(model, newdata = data) # predict lp from the bootstrap model in the original sample
    
    # calculate the apparent performance of the bootstrap model in the bs sample
    app_cstat_model <- roc(outcome~pr_bs,data=bs_samp)
    results[i,1] <- app_cstat_model$auc
    app_citl_model <- glm(outcome ~ offset(lp_bs),family=binomial, data=bs_samp)
    results[i,2] <- summary(app_citl_model)$coefficients[1,1]
    app_cslope_model <- glm(outcome ~ lp_bs,family=binomial(link='logit'), data=bs_samp)
    results[i,3] <- summary(app_cslope_model)$coefficients[2,1]
    
    # calculate the test performance of the bootstrap model in the original sample
    test_cstat_model <- roc(outcome~pr_test,data=data)
    results[i,4] <- test_cstat_model$auc
    test_citl_model <- glm(outcome ~ offset(lp_test),family=binomial, data=data)
    results[i,5] <- summary(test_citl_model)$coefficients[1,1]
    test_cslope_model <- glm(outcome ~ lp_test,family=binomial, data=data)
    results[i,6] <- summary(test_cslope_model)$coefficients[2,1]
    
    print(i)
  }
  results2 <- as.data.frame(results)
  colnames(results2) <- c("app_c_stat" ,"app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope")  
  return(results2)
}
boot_results <- manual_boot(smart,50)

# optimism adjusted statistics
c1$auc - (mean(boot_results$app_c_stat)-mean(boot_results$test_c_stat)) # c-stat
mod_log_1$coef - (mean(boot_results$app_citl)-mean(boot_results$test_citl)) # citl
mod_log_2$coef[2]- (mean(boot_results$app_c_slope)-mean(boot_results$test_c_slope)) # c-slope

#optimism on its own
mean(boot_results$app_c_stat)-mean(boot_results$test_c_stat)
mean(boot_results$app_citl)-mean(boot_results$test_citl)
mean(boot_results$app_c_slope)-mean(boot_results$test_c_slope)
# optimism adjusted statistics
c1$auc - (mean(boot_results$app_c_stat)-mean(boot_results$test_c_stat)) # c-stat
mod_log_1$coef - (mean(boot_results$app_citl)-mean(boot_results$test_citl)) # citl
mod_log_2$coef[2]- (mean(boot_results$app_c_slope)-mean(boot_results$test_c_slope)) # c-slope

#optimism on its own
mean(boot_results$app_c_stat)-mean(boot_results$test_c_stat)
mean(boot_results$app_citl)-mean(boot_results$test_citl)
mean(boot_results$app_c_slope)-mean(boot_results$test_c_slope)
## To get bootstrapped answers for CITL, c-statistic & c-slope it's a more manual approach:
k10 <- qchisq(0.10,1,lower.tail=FALSE)

manual_boot_bw <- function(data,samples){
  results <- matrix(nrow = samples,ncol = 6)
  set.seed(231398)
  for (i in 1:samples) {
    samp_index <- sample(1:nrow(data), nrow(data), rep=TRUE) # create a sampling index vector
    
    bs_samp <- data[samp_index,] # index the orignal dataset using the sampling vector to give the bs sample
    model <- glm(outcome~SEX+age3_spline + SBP + alcohol + CHOLO + BMIO + DIABETES + CARDIAC + SMOKING + AAA,family=binomial, data=bs_samp)  # Fit model to the bootstrap sample
    backward_model <- stepAIC(model,k=k10,direction="backward",trace=FALSE) 
    pr_bs <- predict(backward_model,type="response")  # predict probabilities from the bootstrap model in the bs sample
    lp_bs <- predict(backward_model) # predict lp from the bootstrap model in the bs sample
    
    pr_test <- predict(model,type="response",newdata = data)  # predict probabilities from the bootstrap model in the original sample
    lp_test <- predict(model, newdata = data) # predict lp from the bootstrap model in the original sample
    
    # calculate the apparent performance of the bootstrap model in the bs sample
    app_cstat_model <- roc(outcome~pr_bs,data=bs_samp)
    results[i,1] <- app_cstat_model$auc
    app_citl_model <- glm(outcome ~ offset(lp_bs),family=binomial, data=bs_samp)
    results[i,2] <- summary(app_citl_model)$coefficients[1,1]
    app_cslope_model <- glm(outcome ~ lp_bs,family=binomial(link='logit'), data=bs_samp)
    results[i,3] <- summary(app_cslope_model)$coefficients[2,1]
    
    # calculate the test performance of the bootstrap model in the original sample
    test_cstat_model <- roc(outcome~pr_test,data=data)
    results[i,4] <- test_cstat_model$auc
    test_citl_model <- glm(outcome ~ offset(lp_test),family=binomial, data=data)
    results[i,5] <- summary(test_citl_model)$coefficients[1,1]
    test_cslope_model <- glm(outcome ~ lp_test,family=binomial, data=data)
    results[i,6] <- summary(test_cslope_model)$coefficients[2,1]
    
    print(i)
  }
  results2 <- as.data.frame(results)
  colnames(results2) <- c("app_c_stat","app_citl","app_c_slope","test_c_stat","test_citl","test_c_slope") 
  return(results2)
}
boot_results_bw <- manual_boot_bw(smart,50)

# optimism adjusted statistics
c1$auc - (mean(boot_results_bw$app_c_stat)-mean(boot_results_bw$test_c_stat)) # c-stat
mod_log_1$coef - (mean(boot_results_bw$app_citl)-mean(boot_results_bw$test_citl)) # citl
mod_log_2$coef[2]- (mean(boot_results_bw$app_c_slope)-mean(boot_results_bw$test_c_slope)) # c-slope

#optimism on its own
mean(boot_results_bw$app_c_stat)-mean(boot_results_bw$test_c_stat)
mean(boot_results_bw$app_citl)-mean(boot_results_bw$test_citl)
mean(boot_results_bw$app_c_slope)-mean(boot_results_bw$test_c_slope)


