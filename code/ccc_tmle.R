# case-control weighted TMLE estimator of Rose and van der Laan (AJE)
# this code demonstrates the double-robustness property of TMLE
# estimating risk in case-control studies
# target parameter: marginal risk difference
library(dplyr)

# using large sample so that small sample/random sampling issues are minimized
#set.seed(1232123)
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

# STEP 1: ASSIGN WEIGHTS
# case-control weights of van der Laan
ccdat$wts = ifelse(ccdat$y, q0, (1-q0)/J)

# STEP 2: WEIGHTED OUTCOME REGRESSION
# weighted logistic regression in case-control data: note this is deliberately misspecified to demonstrate doubly-robust property
(m_theta <- glm(y~x, data=ccdat, family=binomial(), weights=wts))
# initial estimate of outcome model parameters (already given by mtheta, but can also use intercept adjustment)
Q_2initm = m_theta  # distribution of Y | X, Z
# initial prediction logit(P_{X,n}(Y=1|X,Z))
ccdat$Q_2pllink = predict(Q_2initm, type="link")

# STEP 3: WEIGHTED EXPOSURE REGRESSION
# weighted estimate of propensity score parameters (unbiased under this correct model)
(m_eta <- glm(x~z, data = ccdat, family=binomial(), weights=wts))

# STEP 4: CLEVER COVARIATE 
# estimating fluctuation parameter using a clever covariate
g1 = predict(m_eta, type="response")
ccdat$H = ccdat$x/g1 - (1-ccdat$x)/(1-g1)  # clever covariate for risk difference
ccdat$H0 = ccdat$x/g1   # clever covariates for risk ratio
ccdat$H1 = (1-ccdat$x)/(1-g1)  # clever covariates for risk ratio

# STEP 5: UPDATE FIT
# case-control weighted maximum likelihood fit for the submodel
(flucmod = glm(y~-1+H, offset=Q_2pllink, data=ccdat, family=binomial(), weights=wts))
eps = coef(flucmod)[1]     # fluctuation parameter
# fluctuation model for risk ratio
(flucmod2 = glm(y~-1+H0+H1, offset=Q_2pllink, data=ccdat, family=binomial(), weights=wts))
eps2 = coef(flucmod2)[1:2]     # fluctuation parameter


# predictions under initial fit
ccdat$Q_2initlink1 = predict(Q_2initm, newdata=mutate(ccdat, x=1), type="link")
ccdat$Q_2initlink0 = predict(Q_2initm, newdata=mutate(ccdat, x=0), type="link")

# updated predictions, risk difference
ccdat$Q_2uplink1   = ccdat$Q_2initlink1 + eps*1/g1 
ccdat$Q_2uplink0   = ccdat$Q_2initlink0 + eps*(-1/(1-g1))

ccdat$Q_2uplink1rr   = ccdat$Q_2initlink1 + eps2[1]*1/g1 
ccdat$Q_2uplink0rr   = ccdat$Q_2initlink0 + eps2[2]*(1/(1-g1))


# STEP 6: PLUG IN ESTIMATORS OF RISK DIFFERENCE 
# "naive" predictions under exposed/unexposed 
ccdat$Q_2initp1 = plogis(ccdat$Q_2initlink1)
ccdat$Q_2initp0 = plogis(ccdat$Q_2initlink0)
# updated predictions under exposed/unexposed 
ccdat$Q_2upp1 = plogis(ccdat$Q_2uplink1)
ccdat$Q_2upp0 = plogis(ccdat$Q_2uplink0)
# updated predictions under exposed/unexposed (risk ratio)
ccdat$Q_2upp1rr = plogis(ccdat$Q_2uplink1rr)
ccdat$Q_2upp0rr = plogis(ccdat$Q_2uplink0rr)

# standardize over weighted distribution of W via case-control weights
# this corresponds to a case-control weighted average over the sample
# final weighted risk difference
(tmle_estimate <- mean((ccdat$Q_2upp1 - ccdat$Q_2upp0)*ccdat$wts/mean(ccdat$wts)))
(tmle_estimate_rr <- mean((ccdat$Q_2upp1rr)*ccdat$wts/mean(ccdat$wts))/mean((ccdat$Q_2upp0rr)*ccdat$wts/mean(ccdat$wts)))

# compare with naive estimate under misspecified model
(naive_estimate <- mean((ccdat$Q_2initp1 - ccdat$Q_2initp0)*ccdat$wts/mean(ccdat$wts)))
(naive_estimaterr <- mean((ccdat$Q_2initp1)*ccdat$wts/mean(ccdat$wts))/mean((ccdat$Q_2initp0)*ccdat$wts/mean(ccdat$wts)))

##### marginal risk difference estimates
# truth, under true model: E_Z(Y|X=1,Z) - E_Z(Y|X=0,Z)
zvec = seq(0,1, length=1000000)
sum(plogis(-3.5 + 1+zvec+1*zvec)*1/length(zvec)) - sum(plogis(-3.5 + 0+zvec+0*zvec)*1/length(zvec))

# True risk ratio
sum(plogis(-3.5 + 1+zvec+1*zvec)*1/length(zvec))/ (sum(plogis(-3.5 + 0+zvec+0*zvec)*1/length(zvec)))


### bootstrap variance estimation (risk difference only)
#1: create a function that takes in a dataset and outputs an estimate

bootest = function(dat){
  dat = select(dat, c(x,y,z))
  # STEP 1: ASSIGN WEIGHTS
  dat$wts = ifelse(dat$y, q0, (1-q0)/J)
  # STEP 2: WEIGHTED OUTCOME REGRESSION
  (m_theta <- glm(y~x+z, data=dat, family=binomial(), weights=wts))
  Q_2initm = m_theta  # distribution of Y | X, Z
  dat$Q_2pllink = predict(Q_2initm, type="link")
  # STEP 3: WEIGHTED EXPOSURE REGRESSION
  (m_eta <- glm(x~z, data = dat, family=binomial(), weights=wts))
  # STEP 4: CLEVER COVARIATE 
  g1 = predict(m_eta, type="response")
  dat$H = dat$x/g1 - (1-dat$x)/(1-g1)  # clever covariate for risk difference
  # STEP 5: UPDATE FIT
  (flucmod = glm(y~-1+H, offset=Q_2pllink, data=dat, family=binomial(), weights=wts))
  eps = coef(flucmod)[1]     # fluctuation parameter
  dat$Q_2uplink1   = predict(Q_2initm, newdata=mutate(dat, x=1), type="link") + eps*1/g1 
  dat$Q_2uplink0   = predict(Q_2initm, newdata=mutate(dat, x=0), type="link") + eps*1/(1-g1)
  # STEP 6: PLUG IN ESTIMATORS OF RISK DIFFERENCE 
  dat$Q_2upp1 = plogis(dat$Q_2uplink1)
  dat$Q_2upp0 = plogis(dat$Q_2uplink0)
  # final weighted risk difference
  (tmle_estimate <- mean((dat$Q_2upp1 - dat$Q_2upp0)*dat$wts/mean(dat$wts)))
}

#2: create function that samples the data based on case status
bootcc <- function(dat){
  ncases = sum(dat$y)
  ncontrols = nrow(dat)- ncases
  # sample cases
  caseindex = which(dat$y ==1)
  bootcaseindex = sample(caseindex, ncases, replace=TRUE)
  
  # sample controls
  controlindex = which(dat$y ==0)
  bootconttrolindex = sample(controlindex, ncontrols, replace=TRUE)
  rbind(
    dat[bootcaseindex,],
    dat[bootconttrolindex,]
  )
}

# create function that samples and runs analysis
run_bootsample <- function(i){
  bootdat = bootcc(ccdat)
  suppressWarnings(bootest(bootdat))
}

#3: run that function over numerous draws from the data
# not actually run because this is time intensive with large samples
#bootests = sapply(1:1000, run_bootsample)
#bootvar = var(bootests)                 # bootstrap variance
#quantile(bootests, c(0.025, 0.975))     # percentile based confidence limits
#tmle_estimate
