# Corrected intercept approach: attributable to multiple authors
# estimating risk in case-control studies
# target parameters: logistic model parameters, marginal risk difference and conditional risk difference
library(dplyr)

# using large sample so that small sample/random sampling issues are minimized
set.seed(1232123)
N = 1000000
# generate syntetic data under a known logistic model
z = runif(N)                             # confounder
x = rbinom(N, 1, plogis(-.5 + z))        # exposure
y = rbinom(N, 1, plogis(-3.5 + x+z+x*z)) # true model for outcome

# case-control sampling probabilistic 2:1 case:control sampling, take all cases
samp = y==1 | as.logical(rbinom(N, 1, 2*mean(y)/mean(1-y)))

# cohort data
codat = data.frame(z,x,y)
# case control data
ccdat = data.frame(z = z[samp], x=x[samp],y=y[samp])

q0 = mean(codat$y)                                    # population risk
q1 = mean(ccdat$y)                                    # study "risk"
# intercept correction: population log-odds

ccdat$offset = -log((q0/(1-q0))/(q1/(1-q1)))

##### logistic model parameters
# unmodified logistic regression in cohort data (all parameters unbiased)
summary((m1 <- glm(y~x+z+x*z, data=codat, family=binomial())))
# unmodified logistic regression in case-control data (all parameters unbiased except intercept)
summary((m2 <- glm(y~x+z+x*z, data=ccdat, family=binomial())))
# offset logistic regression in case-control data (all parameters unbiased)
summary((m_interceptcorrected <- glm(y~x+z+x*z+offset(offset), data=ccdat, family=binomial())))

# difference in distribution of Z across case-control and cohort data
mean(codat$z) # unbiased
mean(ccdat$z) # biased because cases are oversampled, thus over-sampling Z

##### marginal risk difference estimates
# truth, under true model: E_Z(Y|X=1,Z) - E_Z(Y|X=0,Z)
zvec = seq(0,1, length=1000000)
sum(plogis(-3.5 + 1+zvec+1*zvec)*1/length(zvec)) - sum(plogis(-3.5 + 0+zvec+0*zvec)*1/length(zvec))


# unmodified logistic model in cohort data: unbiased conditional risks, unbiased distribution of confounder
round(mean(predict(m1, newdata = mutate(codat, x=1), type="response") - predict(m1, newdata = mutate(codat, x=0), type="response")), 3)
# unmodified logistic model in case-control data: biased conditional risks, biased distribution of confounder
round(mean(predict(m2, newdata = mutate(ccdat, x=1), type="response") - predict(m2, newdata = mutate(ccdat, x=0), type="response")), 3)
# intercept correction model in case-control data: unbiased conditional risks, biased distribution of confounder
round(mean(plogis(predict(m_interceptcorrected, newdata = mutate(ccdat, x=1), type="link") - ccdat$offset) - plogis(predict(m_interceptcorrected, newdata = mutate(ccdat, x=0), type="link") - ccdat$offset)), 3)


#### conditional risk difference estimates (conditional on z=0.5)
# truth, under true model: E_Z(Y|X=1,Z=0.5) - E_Z(Y|X=0,Z=0.5)
zvec = 0.5
round(sum(plogis(-3.5 + 1+zvec+1*zvec)*1/length(zvec)) - sum(plogis(-3.5 + 0+zvec+0*zvec)*1/length(zvec)), 3)

# cannot actually condition on z=0.5 in the data, so we pick a small window around 0.5
# unmodified logistic model in cohort data: unbiased conditional risks (conditional on Z)
round(mean(predict(m1, newdata = mutate(codat, x=1) %>% filter(z>0.45 & z<0.55), type="response") - predict(m1, newdata = mutate(codat, x=0) %>% filter(z>0.45 & z<0.55), type="response")), 3)
# unmodified logistic model in case-control data: biased conditional risks (conditional on Z)
round(mean(predict(m2, newdata = mutate(ccdat, x=1) %>% filter(z>0.45 & z<0.55), type="response") - predict(m2, newdata = mutate(ccdat, x=0) %>% filter(z>0.45 & z<0.55), type="response")), 3)
# intercept correction model in case-control data: unbiased conditional risks (conditional on Z)
round(mean(plogis(predict(m_interceptcorrected, newdata = mutate(ccdat, x=1) %>% filter(z>0.45 & z<0.55), type="link") - filter(ccdat, z>0.45 & z<0.55)$offset) - plogis(predict(m_interceptcorrected, newdata = mutate(ccdat, x=0) %>% filter(z>0.45 & z<0.55), type="link") - filter(ccdat, z>0.45 & z<0.55)$offset)), 3)
