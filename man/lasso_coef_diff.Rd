\name{lasso_coef_diff}
\alias{lasso_coef_diff}
\title{Importance statistics based on Lasso coefficient}
\description{
This function provides importance scores for variables (including the knockoffs) in order to compute statistics. The variable importance is measured by Lasso coefficient.
}
\usage{
lasso_coef_diff(Xaug, y, family = Gaussian())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xaug}{augmented design matrix combining original predictors and knockoff variables}
  \item{y}{response vector, or a survival object with two columns}
  \item{family}{Binomial(), Binomial(link = “logit”, type=”glm”), Gaussian(), Poisson(), CoxPH(), Cindex(), GammaReg(), NBinomial(), Weibull(), Loglog(), Lognormal(), etc. See mboost documentation for details.}
}
\value{
2p vector containing varible importance for both orginal variables and knockoff variables
}
