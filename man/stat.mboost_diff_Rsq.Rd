\name{stat.mboost_diff_Rsq}
\alias{stat.mboost_diff_Rsq}
\title{Importance statistics based on R-square via boosting}
\description{
This function provides importance scores for variables (including the knockoffs) in order to compute statistics. The variable importance is measured by R-squares obtained from boosting.
}
\usage{
stat.mboost_diff_Rsq(Xaug, y, max.mstop = 100, bl = c("bbs", "bols", "btree"), cv.fold = 5, family = Gaussian())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xaug}{augmented design matrix combining original predictors and knockoff variables}
  \item{y}{response vector, or a survival object with two columns}
  \item{max.mstop}{maximum number of boosting iteration}
  \item{bl}{base-learners when fitting models using mboost. 'bols' means linear base-learners, 'bbs' penalized regression splines with a B-spline basis, and 'btree' boosts stumps.}
  \item{cv.fold}{number of folds in cross-validation to choose number of iteration}
  \item{family}{Binomial(), Binomial(link = “logit”, type=”glm”), Gaussian(), Poisson(), CoxPH(), Cindex(), GammaReg(), NBinomial(), Weibull(), Loglog(), Lognormal(), etc. See mboost documentation for details.}
}
\value{
2p vector containing varible importance for both orginal variables and knockoff variables
}
