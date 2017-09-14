
selection <- function(X, y, response.type = c("continuous","binomial","survival"),
                          q=0.10, knockoff.method = c("sdp","asdp"), knockoff.shrink = T,
                          screen = T, screening.num = nrow(X), screening.knot = 10,
                          max.mstop = 100, baselearner = c("bbs", "bols", "btree"), cv.fold = 5, family = NULL,
                          threshold=c('knockoff','knockoff+')) {
  library(MFKnockoffs)
  library(doMC)

  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n | nrow(y) == n)
  response.type = match.arg(response.type)

  ## create knockoffs
  X = scale(X)
  X_k = MFKnockoffs.create.approximate_gaussian(X, method = match.arg(knockoff.method), shrink = knockoff.shrink)
  Xaug = cbind(X, X_k)   # augmented design matrix

  ## screening and statistics
  if(screen){
    SIS = sis(Xaug, y, screening.num, screening.knot, response.type)   # return 0/1 indicating screened results
    Xaug_sis = Xaug[,SIS==1]

    Z = rep(0,2*p)
    Z[SIS==1] = stat.mboost_varimp2(Xaug_sis, y, response.type, max.mstop, match.arg(baselearner), cv.fold, family)
  }
  else Z = stat.mboost_varimp2(Xaug, y, response.type, max.mstop, baselearner, cv.fold, family)

  ## statistics
  orig = 1:p
  W = Z[orig] - Z[orig+p]

  ## run the knockoff filter
  T = MFKnockoffs.threshold(W, q, match.arg(threshold))   # threshold

  selected = (1:p)[W >= T]

}







