\name{selection}
\alias{selection}
\title{Generalized Model-Free Knockoff Filter for Controlled Variable Selection}
\description{
This function is the main entry point for the package. It runs the knockoff procedure which accomodates various covariate distributions and model-free associations using gradient boosting, and ultimately provides the selected variables for the input dataset. Users can specify whether to screen variables before selection and the extent of screening.
}
\usage{
knockoffs.varsel(X, y, family = Gaussian(),
  q=0.10, knockoff.method = c("sdp","asdp"), knockoff.shrink = T,
  screen = T, screening.num = nrow(X), screening.knot = 10,
  max.mstop = 100, baselearner = c("bbs", "bols", "btree"), cv.fold = 5,
  threshold=c('knockoff','knockoff+'))
}
\arguments{
  \item{X}{matrix of predictors}
  \item{y}{response vector, or a survival object with two columns}
  \item{family}{Binomial(), Binomial(link = “logit”, type=”glm”), Gaussian(), Poisson(), CoxPH(), Cindex(), GammaReg(), NBinomial(), Weibull(), Loglog(), Lognormal(), etc. See mboost documentation for details.}
  \item{q}{target FDR (false discovery rate)}
  \item{knockoff.method}{method to construct knockoffs, either SDP or approximate SDP}
  \item{knockoff.shrink}{whether to shrink the estimated covariance matrix (default: FALSE)}
  \item{screen}{whether to screen the variables (default: TRUE)}
  \item{screening.num}{number of variables left after screening (default: sample size)}
  \item{screening.knot}{parameter in screening process}
  \item{max.mstop}{maximum number of boosting iteration}
  \item{baselearner}{base-learners when fitting models using mboost. 'bols' means linear base-learners, 'bbs' penalized regression splines with a B-spline basis, and 'btree' boosts stumps.}
  \item{cv.fold}{number of folds in cross-validation to choose number of iteration}
  \item{threshold}{method to calculate knockoff threshold, either “knockoff” or “knockoff+”}
}
\details{
The default family for continuous response is Gaussian(), Binomial() for binary response, and CoxPH() for survival response.
}
\value{
A vector containing the selected covariate indices
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
