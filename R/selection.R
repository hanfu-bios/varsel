
selection <- function(X, y, family = Gaussian(),
                          q=0.10, knockoff.method = c("sdp","asdp"), knockoff.shrink = T,
                          screen = T, screening.num = nrow(X), screening.knot = 10,
                          max.mstop = 100, baselearner = c("bbs", "bols", "btree"), cv.fold = 5,
                          threshold=c('knockoff','knockoff+')) {
  library(MFKnockoffs)
  library(doMC)

  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n | nrow(y) == n)


  ## create knockoffs
  X = scale(X)
  X_k = MFKnockoffs.create.approximate_gaussian(X, method = match.arg(knockoff.method), shrink = knockoff.shrink)
  Xaug = cbind(X, X_k)   # augmented design matrix

  ## screening and statistics
  if(screen){
    check.dist <- function(x) grepl(x,family@name)
    if (any(unlist(lapply(c("Cox","Weibull","Log Logistic","Lognormal","Gehan","Concordance Probability"),check.dist))))  #survival
        SIS = sis.surv(Xaug, y, screening.num, screening.knot)   # return 0/1 indicating screened results
    else SIS = sis(Xaug, y, screening.num, screening.knot)
    Xaug_sis = Xaug[,SIS==1]

    Z = rep(0,2*p)
    Z[SIS==1] = stat.mboost_varimp2(Xaug_sis, y, max.mstop, match.arg(baselearner), cv.fold, family)
  }
  else Z = stat.mboost_varimp2(Xaug, y, max.mstop, baselearner, cv.fold, family)

  ## statistics
  orig = 1:p
  W = Z[orig] - Z[orig+p]

  ## run the knockoff filter
  T = MFKnockoffs.threshold(W, q, match.arg(threshold))   # threshold

  selected = (1:p)[W >= T]

}







