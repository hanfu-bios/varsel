library(knockoffs.varsel)
rm(list=ls())
set.seed(1234)
n = 100; p = 30; k = 10; amplitude = 1    # k is the number of true signals

# Multivariate normal design matrix X with AR(1) covariance
mu = rep(0,p); rho = 0.3
Sigma = diag(p); for (a in 1:p){for(b in 1:p) Sigma[a,b]=rho^abs(a-b)}
library(MASS)
X.sample <- function() mvrnorm(n, mu=rep(0,p), Sigma=Sigma)

# Exponential
# X.sample <- function() {
#   Xnorm = mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
#   X = matrix(0,n,p)
#   for (j in 1:p) X[,j] = qexp(rank(Xnorm[,j],ties.method = "average")/(n+1))
#   X
# }

# beta's
nonzero = sample(p, k)    # true signal indices
beta.sample <- function(amplitude){
  beta = rep(0,p)
  beta[1:p %in% nonzero] = amplitude * sample(rep(c(1,-1), length.out = k))
  return(beta)
}

fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
power <- function(selected) sum(beta[selected] != 0) / k

X = X.sample()
beta = beta.sample(amplitude)

# gaussian response
y = X^2 %*% beta + rnorm(n)
selected = selection(X,y)
fdp(selected)
power(selected)

# binomial response

y = rbinom(n, 1, exp(X %*% beta)/(1+exp(X %*% beta)))
selected = selection(X,y, family = Binomial())
fdp(selected)
power(selected)

# logistic regression
selected = selection(X,y, family = Binomial(link = "logit"), baselearner = "bols")
fdp(selected)
power(selected)


# survival response 1
U=runif(n, 0,1)
pre_time=-log(U)/(1*exp(X %*% beta))
pre_censoring=runif(n,1,30)
pre_censoring=pre_censoring*(pre_censoring<3)+3*(pre_censoring>=3)
tcens=(pre_censoring<pre_time) # censoring indicator
delta=1-tcens
time=pre_time*(delta==1)+pre_censoring*(delta==0)
library(survival)
y = Surv(time, delta)
selected = selection(X,y, family = CoxPH())
fdp(selected)
power(selected)

# survival response 2
library(survival)
T = rweibull(n, shape=1, scale=.002*exp(X %*% beta))   # true event time
C = rweibull(n, shape=1, scale=.004)   #censoring time
time = pmin(T,C)  #observed time is min of censored and true
event = time==T
y = Surv(time, event)
selected = selection(X,y, family = CoxPH())
fdp(selected)
power(selected)
