# case-control weighted ML estimator of van der Laan vs. model offset approach
# estimating risk in case-control studies
# target parameter: marginal risk difference
library(dplyr)

# using large sample so that small sample/random sampling issues are minimized
set.seed(1232123)
N = 1000000
# generate syntetic data under a known logistic model
z = runif(N)                    # confounder
x = rbinom(N, 1, plogis(-.5 + z))        # exposure
y = rbinom(N, 1, plogis(-3.5 + x+z+x*z)) # true model for outcome

# case-control sampling probabilistic 2:1 case:control sampling, take all cases
samp = y==1 | as.logical(rbinom(N, 1, 2*mean(y)/mean(1-y)))

# cohort data
codat = data.frame(z,x,y)
# case control data
ccdat = data.frame(z = z[samp], x=x[samp],y=y[samp])

(J = sum(1-ccdat$y)/sum(ccdat$y))                     # control:case ratio in case-control data
q0 = mean(codat$y)                                    # population risk

# case-control weights of van der Laan
ccdat$wts = ifelse(ccdat$y, q0, (1-q0)/J)

##### logistic model parameters
# unmodified logistic regression in cohort data (all parameters unbiased)
summary((m1 <- glm(y~x+z+x*z, data=codat, family=binomial())))
# unmodified logistic regression in case-control data (all parameters unbiased except intercept)
summary((m2 <- glm(y~x+z+x*z, data=ccdat, family=binomial())))
# weighted logistic regression in case-control data (all parameters unbiased)
summary((m_weighted <- glm(y~x+z+x*z, data=ccdat, family=binomial(), weights=wts)))


##### marginal risk difference estimates
# truth, under true model: E_Z(Y|X=1,Z) - E_Z(Y|X=0,Z)
zvec = seq(0,1, length=1000)
sum(plogis(-3.5 + 1+zvec+1*zvec)*1/length(zvec)) - sum(plogis(-3.5 + 0+zvec+0*zvec)*1/length(zvec))


# unbiased in unmodified cohort model: unbiased conditional risks, unbiased distribution of confounder
round(mean(predict(m1, newdata = mutate(codat, x=1), type="response") - predict(m1, newdata = mutate(codat, x=0), type="response")), 3)
# biased in unmodified case-control model: biased conditional risks, biased distribution of confounder
round(mean(predict(m2, newdata = mutate(ccdat, x=1), type="response") - predict(m2, newdata = mutate(ccdat, x=0), type="response")), 3)
# unbiased in weighted case-control model: unbiased weighted conditional risks, unbiased weighted distribution of confounder
round(mean(predict(m_weighted, newdata = mutate(ccdat, x=1), type="response")*ccdat$wts/mean(ccdat$wts)) - mean(predict(m_weighted, newdata = mutate(ccdat, x=0), type="response")*ccdat$wts/mean(ccdat$wts)), 3)

