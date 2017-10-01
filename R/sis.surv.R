sis.surv <- function(Xaug, y, N, knot){
  library(splines)
  n = nrow(Xaug)
  GD=NULL
  select_SIS=rep(0,ncol(Xaug))

    delta = y[,2]; time = y[,1];
    delta = delta[order(time)]
    Xaug = Xaug[order(time),]
    time = time[order(time)]
    for (j in 1:(ncol(Xaug))){
      knot_set=quantile(Xaug[,j],prob=seq(1:(knot-4))/(knot-3))
      bs7=bs(Xaug[,j],df=knot, knot=knot_set, intercept=FALSE, degree=3)
      bs8=matrix(bs7, nrow=n) # for build beta_t
      beta=rep(0,knot-1)
      result=ddloglik(n,delta,bs8,beta)
      #L1_summary=cbind(L1_summary,result$L1)
      #L2_summary=cbind(L2_summary,result$L2)
      update=qr.solve(result$L2,result$L1,  tol = 1e-20)

      #temp3<-coxph(Surv(time, delta) ~ bs8)
      #update=summary(temp3)$coefficients[,1]
      update=update/sqrt(sum(update^2))
      GD[j]=sum(result$L1*update)
      #z_spline=cbind(z_spline, bs8)
    }
    select_SIS[rev(order(GD))[1:N]]=1

  select_SIS
}
