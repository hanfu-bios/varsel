\name{create.svm}
\alias{create.svm}
\title{Regression-based knockoffs construction
}
\description{This function constructs knockoffs based on regression, specifically support vector regression, and permutation of errors, given the conditionally exchangeable feature of knockoffs.
}
\usage{
create.svm(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of original predictors}
}
\value{
n by p knockoffs matrix
}
