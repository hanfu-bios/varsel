sis <- function(Xaug, y, N, knot){
  library(splines)
  n = nrow(Xaug)
  GD=NULL
  select_SIS=rep(0,ncol(Xaug))
    y = y - mean(y)
    for (j in 1:(ncol(Xaug))){
      knot_set=quantile(Xaug[,j],prob=seq(1:(knot-4))/(knot-3))
      bs7=bs(Xaug[,j],df=knot, knot=knot_set, intercept=TRUE, degree=3)
      bs8=matrix(bs7, nrow=n) # for build beta_t
      #z_spline=cbind(z_spline, bs8)
      GD[j]=t(t(bs8)%*%y)%*%solve(t(bs8)%*%bs8)%*%(t(bs8)%*%y)
    }
    select_SIS[rev(order(GD))[1:N]]=1

  select_SIS
}
