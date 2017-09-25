\name{mboost.varimp}
\alias{mboost.varimp}
\title{Variable Importance for Boosting}
\description{
This function provides importance scores for variables (including the knockoffs) in order to compute statistics. The variable importance is measured by in-bag risk reduction per boosting step.
}
\usage{
mboost.varimp(Xaug, y, max.mstop = 100, baselearner = c("bbs", "bols", "btree"), cv.fold = 5, family = Gaussian())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xaug}{augmented design matrix combining original predictors and knockoff variables}
  \item{y}{response vector, or a survival object with two columns}
  \item{max.mstop}{maximum number of boosting iteration}
  \item{baselearner}{base-learners when fitting models using mboost. 'bols' means linear base-learners, 'bbs' penalized regression splines with a B-spline basis, and 'btree' boosts stumps.}
  \item{cv.fold}{number of folds in cross-validation to choose number of iteration}
  \item{family}{Binomial(), Binomial(link = “logit”, type=”glm”), Gaussian(), Poisson(), CoxPH(), Cindex(), GammaReg(), NBinomial(), Weibull(), Loglog(), Lognormal(), etc. See mboost documentation for details.}
}
\value{
2p vector containing varible importance for both orginal variables and knockoff variables
}
\references{
  Benjamin Hofner, Andreas Mayr, Nikolay Robinzonov and Matthias Schmid (2014). Model-based Boosting in R: A Hands-on Tutorial Using the R Package mboost. Computational Statistics, 29, 3–35.
  http://dx.doi.org/10.1007/s00180-012-0382-5
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
