\name{sis.surv}
\alias{sis.surv}
\title{Variable Screening before Selection for Survival Response}
\description{When the response is a survival object, this function screens the variables before computing their statistics to reduce computing time, especially when p is large.}
\usage{
sis(Xaug, y, N = nrow(Xaug), knot = 10)
}
\arguments{
  \item{Xaug}{augmented design matrix combining original predictors and knockoff variables}
  \item{y}{response vector, or a survival object with two columns}
  \item{N}{number of variables left after screening (default: sample size)}
  \item{knot}{parameter in screening process}
}
\details{
Sure Independence Screening (SIS) method is used in this function. The screening process would accelerate computing but not necessarily improve power as true signals might be removed by this process. Our experiment indicates that the power would remain almost the same.
}
\value{
a binary 2p by 1 vector for screening results. 1 means selected by the screening and 0 otherwise.
}
\references{
Jianqing Fan and Jinchi Lv (2008) Sure Independence Screening for Ultrahigh Dimensional Feature
Space (with discussion). Journal of Royal Statistical Society B, 70, 849-911.
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
