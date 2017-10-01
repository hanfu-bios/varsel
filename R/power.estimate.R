
power.estimate <- function(n, p, X.dist=c("Gaussian","Binary","Exponential"), X.mu=rep(0,p), X.cov=diag(p),
                           beta = NULL, numTrue = NULL, percentTrue = NULL, amplitude=1,
                           association = c("linear","power","exponential","cosine"), power.degree=2,
                           link = c("identity","logit","survival"), family = NULL,
                           surv.lambdaT=.002, surv.lambdaC=.004, surv.shape=1,
                           nIterations = 10, ...){

  library(MASS)

  stopifnot(length(X.mu)==p)
  stopifnot(all(dim(X.cov)==c(p,p)))
  if (is.null(c(beta, numTrue, percentTrue))) stop("coefficient is not specified")
  if (is.null(numTrue)) numTrue = max(round(p*percentTrue), sum(beta==0))

  norm2exp <- function(x) qexp(rank(x,ties.method = "average")/(length(x)+1))
  X.sample <- function(dist){
    X = mvrnorm(n, mu=X.mu, Sigma=X.cov)
    switch(dist,
           "Gaussian" = X,
           "Binary" = matrix(as.numeric(X>0),n,p),
           "Exponential" = apply(X, 2, norm2exp))
  }

  beta.sample <- function(amplitude){

    nonzero = sample(p, numTrue)
    beta = rep(0,p)
    beta[1:p %in% nonzero] = amplitude * sample(c(1,-1), replace = T)
    return(beta)
  }

  y.sample <- function(X, beta, assoc, link, lambdaT, lambdaC, shape) {
    linear.predictor = switch(assoc,
                              "linear" = X %*% beta,
                              "power" = X^power.degree %*% beta,
                              "exponential" = exp(X) %*% beta,
                              "cosine" = cos(X) %*% beta)
    if (link == "survival"){
      library(survival)
      T = rweibull(n, shape=shape, scale=lambdaT*exp(linear.predictor))   # true event time
      C = rweibull(n, shape=shape, scale=lambdaC)   #censoring time
      time = pmin(T,C)  #observed time is min of censored and true
      event = time==T   # set to 1 if event is observed
    }
    switch(link,
           "identity" = linear.predictor,
           "logit" = rbinom(n, 1, exp(linear.predictor)/(1+exp(linear.predictor))),
           "survival" = Surv(time, event))
  }
  fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  power <- function(selected) sum(beta[selected] != 0) / numTrue

  power_list = NULL
  fdp_list = NULL


  for (iter in 1:nIterations){
    X = X.sample(match.arg(X.dist))
    if (is.null(beta)) beta = beta.sample(amplitude)
    y = y.sample(X, beta, match.arg(association), match.arg(link), surv.lambdaT, surv.lambdaC, surv.shape)

    if(is.null(family)) family = switch(match.arg(link),
                                       "identity" = Gaussian(),
                                       "logit" = Binomial(),
                                       "survival" = CoxPH())
    selected = selection(X, y, family = family, ...)
    power_list[iter] = power(selected)
    fdp_list[iter] = fdp(selected)
  }
  list(power.mean = mean(power_list),
       power.all = power_list,
       power.sd = sd(power_list),
       fdp = mean(fdp_list)
  )
}
