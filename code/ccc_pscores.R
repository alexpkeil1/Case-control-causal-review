# NOTE this file is a placeholder and is not complete
# propsensity score based methods for case-control studies
# estimating risk in case-control studies
# target parameters: population level odds ratios, marginal risk difference, conditional causal effects
library(dplyr)

# using large sample so that small sample/random sampling issues are minimized
set.seed(1232123)
N = 1000000
# generate syntetic data under a known logistic model
z = rbinom(N, 1, 0.5)                    # confounder
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

# case-control weights
ccdat$wts = ifelse(ccdat$y, q0, (1-q0)/J)
ccdat$offset = -log(q0/(1-q0))

#####  propensity score models
# 1: fit to controls (invoking rare disease assumption)
(psmod_controls <- glm(x~z, data = filter(ccdat, y==0), family=binomial()))
# 2: fit to weighted sample
(psmod_weighted <- glm(x~z, data = ccdat, family=binomial(), weights=wts))

#####  propensity scores models
ps_controls <- predict(psmod_controls, type="response")
ps_weighted <- predict(psmod_weighted, type="response")

##### logistic model parameters
# TODO: complete