\name{power.estimate}
\alias{power.estimate}
\title{Power Estimation by Generalized Model-Free Knockoffs Filter}
\description{
This function could be used for estimating power and FDP using the knockoffs filter prior to data collection. Once the user inputs the dimensions of data (sample size and number of covariates), and certain expectation for data structure and association type, this function could simulate data for multiple times and ultimately give an expected value for power and FDR.
}
\usage{
power.estimate(n, p, X.dist=c("Gaussian","Binary","Exponential"), X.mu=rep(0,p), X.cov=diag(p),
                           beta = NULL, numTrue = NULL, percentTrue = NULL, amplitude=1,
                           association = c("linear","power","exponential","cosine"), power.degree=2,
                           link = c("identity","logit","survival"), family = NULL,
                           surv.lambdaT=.002, surv.lambdaC=.004, surv.shape=1,
                           nIterations = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{p}{number of covariates, including null variables}
  \item{X.dist}{distribution of design matrix. Either "Gaussian", "Binary" or "Exponential"}
  \item{X.mu}{expected values for X, a vector of length p (default: zero vector of length p)}
  \item{X.cov}{variance-covariance matrix (p by p) for X (default: identity)}
  \item{beta}{coefficients for p variables if known, a vector of length p}
  \item{numTrue}{number of true signals among p variables}
  \item{percentTrue}{percentage of true signals among p variables}
  \item{amplitude}{signal amplitude}
  \item{association}{association between predictors and response (on the scale of linear predictors). The linear predictor will be X*beta when the input argument is "linear", X^[some power]*beta when "power", exp(X)*beta when "exponential", and cos(X)*beta when "cosine".}
  \item{power.degree}{power degree when the "power" association is selected (default: 2)}
  \item{link}{link function between linear predictor and the response. "identity" for identity link and "logit" for logit link. If "survival" is selected, then survival response will be generated using the hazard function in Cox model.}
  \item{family}{Binomial(), Binomial(link = “logit”, type=”glm”), Gaussian(), Poisson(), CoxPH(), Cindex(), GammaReg(), NBinomial(), Weibull(), Loglog(), Lognormal(), etc. See mboost documentation for details.}
  \item{surv.lambdaT}{baseline hazard in survival response, default: 0.002}
  \item{surv.lambdaC}{hazard of censoring in survival response, default: 0.004}
  \item{surv.shape}{shape parameter of weibull distribution, default: 1}
  \item{nIterations}{number of runs to get the means / distributions of estimated power and FDR}
  \item{...}{further arguments passed to function selection}
}
\details{
At least one of the three arguments, beta, numTrue, and percentTrue, must be specified, or, an error would appear.
For now, the signal amplitude is set to be identical for all the true signals. Generalizations could be made in the future.
}
\value{
A list containing expected value of power, a list of power values from all experiments, standard deviation of power, mean value of FDR achieved (expected to be around the target value)
}
\references{
  Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
  arXiv:1610.02351 (2016).
  \href{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}

  Barber and Candes,
  Controlling the false discovery rate via knockoffs.
  Ann. Statist. 43 (2015), no. 5, 2055--2085.
  \href{https://projecteuclid.org/euclid.aos/1438606853}{https://projecteuclid.org/euclid.aos/1438606853}

  Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid (2014). Model-based Boosting in R: A Hands-on Tutorial Using the R Package mboost. Computational Statistics, 29, 3–35.
  \href{http://dx.doi.org/10.1007/s00180-012-0382-5}{http://dx.doi.org/10.1007/s00180-012-0382-5}
  Available as vignette via: vignette(package = "mboost", "mboost_tutorial")
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
