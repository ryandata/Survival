# Survival Analysis
# Sample Code
# Ryan Womack (libguides.rutgers.edu/data_R) or (github.com/ryandata/Survival)
# Fall 2014 (code rev. 2016-12-23)


#####################
# Section 1
# installing packages
#####################

install.packages("car", dependencies=TRUE)
install.packages("survival", dependencies=TRUE)
install.packages("flexsurv", dependencies=TRUE)
install.packages("KMsurv", dependencies=TRUE)
install.packages("e1071", dependencies=TRUE)
install.packages("rms", dependencies=TRUE)


#####################
# Section 2
# loading and exploring data
#####################

gbcs <- read.csv("https://ryanwomack.com/data/gbcs.csv")
attach(gbcs)
names(gbcs)
dim(gbcs)
summary(gbcs)
hist(age)
plot(density(age))
table(menopause)
table(hormone)
hist(prog_recp)
plot(density(prog_recp), main="Progesterone receptors, density plot")

#Fig 1, prog_recp
par(mfrow=c(1,2))
plot(prog_recp, main="Progesterone receptors")
plot(log(prog_recp), main="log of Progesterone receptors")
par(mfrow=c(1,1))

hist(estrg_recp)
plot(density(estrg_recp))
table(censrec)
hist(rectime)
plot(density(rectime))

#plots
plot(rectime~id)
plot(rectime~age)
plot(rectime~menopause)
plot(rectime~hormone)
plot(rectime~prog_recp)
plot(rectime~estrg_recp)
plot(rectime~censrec)

#correlation matrix
cor(gbcs)

#fancy scatterplot
library(car)
scatterplotMatrix(gbcs)

#####################
# Section 3
# Survival object and survfit
#####################

#reuseable survival objects
library(survival)
recsurv<-Surv(rectime,censrec)


#survfit fits survival curves with various methods
#Kaplan-Meier is most common, so use that if in doubt

fit_KM <- survfit(recsurv~1,type="kaplan-meier", conf.type="log-log")
fit_FH <- survfit(recsurv~1,type="fleming-harrington", conf.type="log-log")
fit_FH2 <- survfit(recsurv~1, type="fh2", conf.type="log-log")

#plot survival
plot(fit_KM, main="survival function for rectime (K-M estimate)", xlab="days", ylab="p")
plot(fit_FH)
plot(fit_FH2)

#print restricted means
print(fit_KM,print.rmean=TRUE)
print(fit_FH,print.rmean=TRUE)
print(fit_FH2,print.rmean=TRUE)

#can plot cumulative hazard, and other mods as 
plot(fit_KM, fun="cumhaz")

#cumulative events (f(y)=1-y)
plot(fit_KM, fun="event")

#log(-log(y))
plot(fit_KM, fun="cloglog")

#log t
plot(fit_KM, fun="log")

#to generate figure for report
par(mfrow=c(1,2))
plot(fit_KM, main="survival function for rectime (K-M estimate)", xlab="days", ylab="p")
plot(fit_KM, fun="cumhaz", main="cumulative hazard function for rectime (K-M estimate)", xlab="days", ylab="p")
par(mfrow=c(1,1))

#survfits to illustrate impact of variables
leg.txt<-c("0", "1")
fit <- survfit(recsurv~as.numeric(age>median(age)))
plot(fit, col=c(2,4))
legend("topright",leg.txt,col=c(2,4),lty=1)

fit <- survfit(recsurv~menopause)
plot(fit, col=c(2,4))
legend("topright",leg.txt,col=c(2,4),lty=1)

fit <- survfit(recsurv~hormone)
plot(fit, col=c(2,4))
legend("topright",leg.txt,col=c(2,4),lty=1)

fit <- survfit(recsurv~as.numeric(prog_recp>median(prog_recp)))
plot(fit, col=c(2,4))
legend("topright",leg.txt,col=c(2,4),lty=1)

fit <- survfit(recsurv~as.numeric(estrg_recp>median(estrg_recp)))
plot(fit, col=c(2,4))
legend("topright",leg.txt,col=c(2,4),lty=1)


#####################
# Section 4
# empirical cdf
#####################

#categorical variables
progtest  <- as.numeric(prog_recp>median(prog_recp))
estrgtest <- as.numeric(estrg_recp>median(estrg_recp))
agetest   <- as.numeric(age>median(age))
                                            
#ecdf plots empirical cdf
plot(ecdf(age[menopause==1]), verticals=TRUE, pch=46, col=2)
lines(ecdf(age[menopause==2]), verticals=TRUE, pch=46, col=4)

plot(ecdf(age[hormone==1]), verticals=TRUE, pch=46, col=2)
lines(ecdf(age[hormone==2]), verticals=TRUE, pch=46, col=4)

plot(ecdf(age[prog_recp>median(prog_recp)]), verticals=TRUE, pch=46, col=2)
lines(ecdf(age[prog_recp<=median(prog_recp)]), verticals=TRUE, pch=46, col=4)

plot(ecdf(age[estrg_recp>median(estrg_recp)]), verticals=TRUE, pch=46, col=2)
lines(ecdf(age[estrg_recp<=median(estrg_recp)]), verticals=TRUE, pch=46, col=4)


#####################
# Section 5
# survreg to check distributions
#####################

#survreg, 1 variable
#here we are checking which distributional assumption provides
#the best fit
fit_exp<-survreg(recsurv~1, dist="exponential")
fit_weibull<-survreg(recsurv~1, dist="weibull")
fit_gauss<-survreg(recsurv~1, dist="gaussian")
fit_logistic<-survreg(recsurv~1, dist="logistic")
fit_lognormal<-survreg(recsurv~1, dist="lognormal")
fit_loglogistic<-survreg(recsurv~1, dist="loglogistic")
summary(fit_exp)
summary(fit_weibull)
summary(fit_gauss)
summary(fit_logistic)
summary(fit_lognormal)
summary(fit_loglogistic)

#try with package "flexsurv"
#flexsurv provides access to additional distributions

library(flexsurv)
fit_exp<-flexsurvreg(recsurv~1, dist="exp")
fit_weibull<-flexsurvreg(recsurv~1, dist="weibull")
fit_gamma<-flexsurvreg(recsurv~1, dist="gamma")
fit_gengamma<-flexsurvreg(recsurv~1, dist="gengamma")
fit_genf<-flexsurvreg(recsurv~1, dist="genf")
fit_lognormal<-flexsurvreg(recsurv~1, dist="lnorm")
fit_gompertz<-flexsurvreg(recsurv~1, dist="gompertz")
fit_exp
fit_weibull
fit_gamma
fit_gengamma
fit_genf
fit_lognormal
fit_gompertz
plot(fit_exp)
plot(fit_weibull)
plot(fit_gamma)
plot(fit_gengamma)
plot(fit_genf)
plot(fit_lognormal)
plot(fit_gompertz)

#log likelihood test
fit_exp$loglik
fit_weibull$loglik
fit_gamma$loglik
fit_gengamma$loglik
fit_lognormal$loglik
fit_gompertz$loglik

qchisq(.95,2)
qchisq(.95,1)

#test gengamma vs. lognormal
f_lognormal<-2*(fit_gengamma$loglik - fit_lognormal$loglik) #2 df
f_lognormal


#AIC is reported by flexsurv
fit_exp$AIC
fit_weibull$AIC
fit_gamma$AIC
fit_gengamma$AIC
fit_lognormal$AIC
fit_gompertz$AIC

#####################
# Section 6
# more exploration
#####################

#KMsurv package can generate life table
#not really needed in this analysis
library(KMsurv)
lifetab()

#log transforms
logprog<-log(prog_recp+.1)
logestrg<-log(estrg_recp+.1)
logrectime<-log(rectime)
logrecsurv<-Surv(logrectime,censrec)

#work with logrec a bit
fit_exp<-flexsurvreg(logrecsurv~1, dist="exp")
fit_weibull<-flexsurvreg(logrecsurv~1, dist="weibull")
fit_gamma<-flexsurvreg(logrecsurv~1, dist="gamma")
fit_gengamma<-flexsurvreg(logrecsurv~1, dist="gengamma")
fit_genf<-flexsurvreg(logrecsurv~1, dist="genf")
fit_lognormal<-flexsurvreg(logrecsurv~1, dist="lnorm")
fit_gompertz<-flexsurvreg(logrecsurv~1, dist="gompertz")
fit_exp
fit_weibull
fit_gamma
fit_gengamma
fit_genf
fit_lognormal
fit_gompertz
plot(fit_exp)
plot(fit_weibull)
plot(fit_gamma)
plot(fit_gengamma)
plot(fit_genf)
plot(fit_lognormal)
plot(fit_gompertz)

#these options work too
plot(fit_gamma, type="hazard")
plot(fit_gamma, type="cumhaz")
plot(fit_gamma, type="survival")

#log likelihood test
fit_exp$loglik
fit_weibull$loglik
fit_gamma$loglik
fit_gengamma$loglik
fit_genf$loglik
fit_lognormal$loglik
fit_gompertz$loglik

qchisq(.95,2)
qchisq(.95,1)

#test genf vs. gengamma
f_gengamma<-2*(fit_genf$loglik-fit_gengamma$loglik) #2 df
f_gengamma

#test genf vs. lognormal
f_lognormal<-2*(fit_genf$loglik - fit_lognormal$loglik) #2 df
f_lognormal

#lognormal is good enough here

#AIC
fit_exp$AIC
fit_weibull$AIC
fit_gamma$AIC
fit_gengamma$AIC
fit_genf$AIC
fit_lognormal$AIC
fit_gompertz$AIC

#confirms choice of lognormal

#create graphic

par(mfrow=c(1,2))
fit_gengamma<-flexsurvreg(recsurv~1, dist="gengamma")
plot(fit_gengamma, main="rectime fit to generalized gamma", xlab="days", ylab="p")

fit_lognormal<-flexsurvreg(logrecsurv~1, dist="lnorm")
plot(fit_lognormal, main="log(rectime) fit to lognormal", xlab="days", ylab="p")
par(mfrow=c(1,1))


#####################
# Section 7
# probability plots with e1071
#####################
#probability plots to test for distribution of response variable
library(e1071)
probplot(rectime)
probplot(rectime, "qunif") #best fit
probplot(rectime, "qexp")
probplot(rectime, "qnorm")
probplot(rectime, "qweibull", shape=1)
probplot(rectime, "qlnorm")
probplot(rectime, "qgamma", shape=1)

#check subgroups, strata

#try different "type" options in coxph


#####################
# Section 8
# CI and diffplots with rms
#####################

#survplot from rms package
#npsurv is equivalent to survreg
#shows CIs easily
library(rms)
fit <- npsurv(recsurv~hormone)
survplot(fit)

#survdiffplot
#shows significant CI
survdiffplot(fit)

#survdiff - numerical test of whether variable influences survival
#rho=0 logrank (Mantel-Haenszel test
#rho=1 Peto & Peto modification of the Gehan-Wilcoxon test

survdiff(recsurv ~ hormone, rho=0)
survdiff(recsurv ~ menopause, rho=0)
survdiff(recsurv ~ agetest, rho=0)
survdiff(recsurv ~ progtest, rho=0)
survdiff(recsurv ~ estrgtest, rho=0)

survdiff(recsurv ~ hormone, rho=1)
survdiff(recsurv ~ menopause, rho=1)
survdiff(recsurv ~ agetest, rho=1)
survdiff(recsurv ~ progtest, rho=1)
survdiff(recsurv ~ estrgtest, rho=1)


#####################
# Section 9
# Cox Proportional Hazards model with coxph and testing fit
#####################

#coxph 

# check each variable
fit <- coxph(recsurv~age)
fit
cox.zph(fit)
plot(cox.zph(fit))
res_martingale<-residuals(fit, type="martingale")
scatter.smooth(age,res_martingale)

fit <- coxph(recsurv~menopause)
fit
cox.zph(fit)
plot(cox.zph(fit))
res_martingale<-residuals(fit, type="martingale")
scatter.smooth(menopause,res_martingale)

fit <- coxph(recsurv~hormone)
fit
cox.zph(fit)
plot(cox.zph(fit))
res_martingale<-residuals(fit, type="martingale")
scatter.smooth(hormone,res_martingale)

fit <- coxph(recsurv~prog_recp)
fit
cox.zph(fit)
plot(cox.zph(fit))
res_martingale<-residuals(fit, type="martingale")
scatter.smooth(prog_recp,res_martingale)

fit <- coxph(recsurv~estrg_recp)
fit
cox.zph(fit)
plot(cox.zph(fit))
res_martingale<-residuals(fit, type="martingale")
scatter.smooth(estrg_recp,res_martingale)

#exact results are not really different than approximate methods
#those methods will allow us to use Schoenfeld residuals
#to check fit with cox.ph
#looking for non-straight line, indicating violation of proportional hazards assumption

#dichotomize estrg_recp 
#still a bad violation, must stratify

fit <- coxph(recsurv ~ estrgtest)
fit
cox.zph(fit)
plot(cox.zph(fit))

#inverse does not converge
fit <- coxph(recsurv ~ I(estrgtest^(-1)))
fit
cox.zph(fit)
plot(cox.zph(fit))

#####################
# Section 10
# residuals and model fitting
#####################

#can also access residuals of other types with
mres<-residuals(fit, type="martingale")
dfb<-residuals(fit, type="dfbetas")
scatter.smooth(rectime,mres)

#get baseline survival
survfit(fit)

#check collinearity
cor(logprog,logestrg)
cor(gbcs)
cor(prog_recp,estrg_recp)

#do stepAIC to select model
#from MASS package, has stepwisemethod for coxph
#stepAIC also accepts forward and backward options for direction

#no interaction, no stratification
library(MASS)
coxbasemodel<-coxph(recsurv~prog_recp)
no_I_no_S<-stepAIC(coxbasemodel, direction="both", scope=list(lower=.~1,upper=.~prog_recp+hormone+age))


#no interaction, stratification by estrg
no_I_estrg<-coxph(recsurv~prog_recp+hormone+age+strata(estrgtest))


#it is very clear that the stratification makes a big difference in log-likelihood
no_I_no_S$loglik
no_I_estrg$loglik

#interaction, no stratification
I_no_S<-stepAIC(coxbasemodel, direction="both", scope=list(lower=.~1,upper=.~hormone+prog_recp+estrg_recp+hormone:prog_recp+hormone:estrg_recp+prog_recp:estrg_recp))


extractAIC(no_I_estrg)
extractAIC(no_I_no_S)
extractAIC(I_no_S)

#check residuals
fit <- no_I_estrg
fit
cox.zph(fit)
plot(cox.zph(fit))

#####################
# Section 11
# more residuals, Cox-Snell
#####################

#can also access residuals of other types with
res_martingale<-residuals(fit, type="martingale")
res_dfbetas<-residuals(fit, type="dfbetas")
res_score<-residuals(fit, type="score")
res_deviance<-residuals(fit, type="deviance")
res_schoenfeld<-residuals(fit, type="schoenfeld")
res_dfbeta<-residuals(fit, type="dfbeta")
res_scaledsch<-residuals(fit, type="scaledsch")
res_partial<-residuals(fit, type="partial")

scatter.smooth(logrectime,res_martingale)
scatter.smooth(logrectime,res_deviance)

# and so on

#cox-snell residuals
res_cox_snell=censrec-res_martingale

fit_cs=survfit(Surv(res_cox_snell,censrec)~1)
Htilde=cumsum(fit_cs$n.event/fit_cs$n.risk)
plot(fit_cs$time,Htilde,type='s',col='blue')
abline(0,1,col='red',lty=2)

#####################
# Section 12
# AFT model using survreg
#####################

#AFT model uses survreg function
#using lognormal distribution according to earlier fit
#use survreg for stepwise
#then flexsurvreg to generate graphs

fit_lognormal<- survreg(logrecsurv~1, dist="lognormal")
aftmodel <- survreg(logrecsurv~hormone+prog_recp+estrg_recp, dist="lognormal")
aftmodel


flexAFT<-flexsurvreg(logrecsurv~hormone+prog_recp, dist="lnorm")
plot(flexAFT, type="survival")
par(mfrow=c(1,2))
plot(flexAFT, type="cumhaz", sub="cumulative hazard for lognormal AFT")
plot(flexAFT, type="hazard")
par(mfrow=c(1,1))
