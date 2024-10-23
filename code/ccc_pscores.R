# NOTE this file is a placeholder and is not complete
# propsensity score based methods for case-control studies
# estimating risk in case-control studies
# target parameters: population level odds ratios, marginal risk difference, conditional causal effects
library(dplyr)
library(survival)
library(modi) # weighted.quantile function
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

# inverse sampling weights (can also do weights of van der Laan)
ccdat$sampwts = ifelse(ccdat$y, 1, sum(1-codat$y)/sum(1-ccdat$y))

#####  propensity score models (except method B using a subcohort)
# Mansson A: reference: ps model fit to cohort data
(psmod_cohort <- glm(x~z, data = codat, family=binomial()))
# Mansson C: fit to weighted sample 
(psmod_weighted <- glm(x~z, data = ccdat, family=binomial(), weights=sampwts))
# Mansson D: fit to controls (invoking rare disease assumption)
(psmod_controls <- glm(x~z, data = filter(ccdat, y==0), family=binomial()))
# Mansson E: fit to case-control data as sample
(psmod_casecon <- glm(x~z, data = ccdat, family=binomial()))
# Mansson F: including case status and using predictions
ccdat$fakecase = ccdat$y
(psmod_modeledcontrol <- glm(x~z+fakecase, data = ccdat, family=binomial()))
# Zhu: zhu's approach, fit among cases and use matching
(psmod_zhu <- glm(x~z, data = filter(ccdat, y==1), family=binomial()))

#####  propensity scores
ps_cohort <- predict(psmod_cohort, type="response")
ps_weighted <- predict(psmod_weighted, type="response")
ps_controls <- predict(psmod_controls, newdata = ccdat, type="response")
ps_casecon <- predict(psmod_casecon, newdata = ccdat, type="response")
ps_modeledcontrol <- predict(psmod_modeledcontrol, newdata = mutate(ccdat, fakecase=0), type="response")
ps_zhu <- predict(psmod_zhu,newdata = ccdat, type="response")

sum(is.na(ps_zhu))

##### inverse probability weights (Mansson methods only)
ipw_cohort = codat$x/ps_cohort + (1-codat$x)/(1-ps_cohort)
ipw_weighted = ccdat$x/ps_weighted + (1-ccdat$x)/(1-ps_weighted)
ipw_controls = ccdat$x/ps_controls + (1-ccdat$x)/(1-ps_controls)
ipw_casecon = ccdat$x/ps_casecon + (1-ccdat$x)/(1-ps_casecon)
ipw_modeledcontrol = ccdat$x/ps_modeledcontrol + (1-ccdat$x)/(1-ps_modeledcontrol)

##### logistic model parameters

### stratification of the propensity score in the quintiles and get summary over strata
#Note: Mansson used Mantel Haenszel OR, but using conditional logit here, which allows weights
# Mansson A: reference: ps model fit to cohort data
codat$ps_cohort_qs = cut(ps_cohort, quantile(ps_cohort, seq(0,1,by=0.2))) # quantile(ps_cohort, seq(0.2,1,by=0.2))
survival::clogit(y~x + strata(ps_cohort_qs), data = codat, method="efron")


# Mansson C: fit to weighted sample 
weighted_qs = c(min(ps_weighted), sapply(seq(0.2,0.8,by=0.2), function(x) modi::weighted.quantile(ps_weighted, ccdat$sampwts, x)), max(ps_weighted))
ccdat$ps_weighted_qs = cut(ps_weighted, weighted_qs, include.lowest = TRUE) # quantile(ps_weighted, seq(0,1,by=0.2))
survival::clogit(y~x + strata(ps_weighted_qs), data = ccdat, method="efron", weights=sampwts)

# Mansson D: fit to controls (invoking rare disease assumption)
ccdat$ps_controls_qs = cut(ps_controls, quantile(ps_controls, seq(0,1,by=0.2)), include.lowest = TRUE) # quantile(ps_cohort, seq(0.2,1,by=0.2))
survival::clogit(y~x + strata(ps_controls_qs), data = ccdat, method="efron")

# Mansson E: fit to case-control data as sample
ccdat$ps_casecon_qs = cut(ps_casecon, quantile(ps_casecon, seq(0,1,by=0.2)), include.lowest = TRUE) # quantile(ps_cohort, seq(0.2,1,by=0.2))
survival::clogit(y~x + strata(ps_casecon_qs), data = ccdat, method="efron")

# Mansson F: including case status and using predictions
ccdat$ps_modeledcontrol_qs = cut(ps_modeledcontrol, quantile(ps_modeledcontrol, seq(0,1,by=0.2))) # quantile(ps_cohort, seq(0.2,1,by=0.2))
survival::clogit(y~x + strata(ps_modeledcontrol_qs), data = ccdat, method="efron")

# Zhu: matching up to 10 controls with each case based on propensity score among cases
zhu_qs = quantile(ps_zhu, seq(0,1,.1))
ccdat$ps_zhu_qs = as.numeric(cut(ps_zhu, zhu_qs, include.lowest = TRUE)) # quantile(ps_weighted, seq(0,1,by=0.2))

# function to create matched sets (create variable with matched set indicator)
unimatch <- function(casestatus, matchvar, caliper=0.01, maxcontrols=10){
  if(length(unique(casestatus)) != 2) stop("casestatus only allows 2 unique values and no missingness")
  if(length(casestatus) != length(matchvar)) stop("casestatus and matchvar should be same length")
  matchid = numeric(length(casestatus))*NA
  ncases = sum(casestatus == max(casestatus))
  caseidx = sample(which(casestatus == max(casestatus)), ncases)
  controlidx = sample(which(casestatus == min(casestatus)))
  controlpool = controlidx
  for(i in 1:ncases){
    matchid[i] = i
    refps = matchvar[caseidx[i]]
    nmatches = 0
    endrec = FALSE
    #while((nmatches<maxcontrols) && !endrec){
    potcontrols = which((matchvar[controlpool] < refps+caliper) & (matchvar[controlpool] > refps-caliper))
    if(length(potcontrols) > maxcontrols) potcontrols = potcontrols[1:maxcontrols]
    idxps = controlpool[potcontrols]
    matchid[which(controlidx %in% idxps)] = i
    controlpool = setdiff(controlpool, controlpool[potcontrols])
    #}
  }
  matchid
}

ccidx = sample(1:dim(ccdat)[1], 100000)
ccdatsub = ccdat[ccidx, ]# subset used because matching is slow
ccdatsub$matchid = unimatch(ccdatsub$y, ps_zhu[ccidx], maxcontrols=10, caliper=0.5) 
sum(table(ccdatsub$matchid))
ccdatsub = ccdatsub[!is.na(ccdatsub$matchid),]
survival::clogit(y~x + strata(matchid), data = ccdatsub, method="efron")



### include propensity score as adjustment variable in logistic regression
# Mansson A: reference: ps model fit to cohort data
summary(glm(y~x+ps_cohort, family=binomial(), data=codat))

# Mansson C: fit to weighted sample 
summary(glm(y~x+ps_weighted, family=binomial(), data=ccdat, weights=ipw_weighted))

# Mansson D: fit to controls (invoking rare disease assumption)
summary(glm(y~x+ps_controls, family=binomial(), data=ccdat))

# Mansson E: fit to case-control data as sample
summary(glm(y~x+ps_casecon, family=binomial(), data=ccdat))

# Mansson F: including case status and using predictions
glm(y~x+ps_modeledcontrol, family=binomial(), data=ccdat)

### marginal structural models
# Mansson A: reference: ps model fit to cohort data
summary(glm(y~x, family=binomial(), data=codat, weights=ipw_cohort))

# Mansson C: fit to weighted sample 
summary(glm(y~x, family=binomial(), data=ccdat, weights=ipw_weighted))

# Mansson D: fit to controls (invoking rare disease assumption)
summary(glm(y~x, family=binomial(), data=ccdat, weights=ipw_controls))

# Mansson E: fit to case-control data as sample
summary(glm(y~x, family=binomial(), data=ccdat, weights=ipw_casecon))

# Mansson F: including case status and using predictions
summary(glm(y~x, family=binomial(), data=ccdat, weights=ipw_modeledcontrol))

