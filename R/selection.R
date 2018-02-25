
selection <- function(X, y, X_k, family = Gaussian(),
                          q=0.10, knockoff.method = c("sdp","asdp","svm"), knockoff.shrink = T,
                          stat = c("RRB", "LCD", "DRS"),
                          screen = T, screening.num = nrow(X), screening.knot = 10,
                          max.mstop = 100, baselearner = c("bbs", "bols", "btree"), cv.fold = 5,
                          threshold=c('knockoff','knockoff+')) {

  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n | nrow(y) == n)


  ## create knockoffs
  X = scale(X)
  if (match.arg(knockoff.method)=="svm")    X_k = create.svm(X)
  else X_k = create.second_order(X, method = match.arg(knockoff.method), shrink = knockoff.shrink)
  Xaug = cbind(X, X_k)   # augmented design matrix
  Xaug = matrix(Xaug, ncol = ncol(Xaug), dimnames = NULL)

  ## screening
  if(n>p*2 & screen) {screen = F; warning("Screening not performed for low dimensinality (p<n/2)")}
  if(screen){
    check.dist <- function(x) grepl(x,family@name)
    if (any(unlist(lapply(c("Cox","Weibull","Log Logistic","Lognormal","Gehan","Concordance Probability"),check.dist))))  #survival
        SIS = sis.surv(Xaug, y, screening.num, screening.knot)   # return 0/1 indicating screened results
    else SIS = sis(Xaug, y, screening.num, screening.knot)
    Xaug_sis = Xaug[,SIS==1]
  }

  ## statistics
  Z = rep(0,2*p)

  if (match.arg(stat) == "RRB"){
    if(screen) Z[SIS==1] = stat.mboost_varimp2(Xaug_sis, y, max.mstop, match.arg(baselearner), cv.fold, family)
    else Z = stat.mboost_varimp2(Xaug, y, max.mstop, match.arg(baselearner), cv.fold, family)
  }

  else if(match.arg(stat) == "DRS"){
    if(screen) Z[SIS==1] = stat.mboost_diff_Rsq(Xaug_sis, y, max.mstop, match.arg(baselearner), cv.fold, family)
    else Z = stat.mboost_diff_Rsq(Xaug, y, max.mstop, match.arg(baselearner), cv.fold, family)
  }

  else if(match.arg(stat) == "LCD"){
    if(screen) Z[SIS==1] = lasso_coef_diff(Xaug_sis, y, family)
    else Z = lasso_coef_diff(Xaug, y, family)
  }

  ## statistics
  orig = 1:p
  W = Z[orig] - Z[orig+p]

  ## run the knockoff filter
  offset = ifelse(match.arg(threshold)=="knockoff", 0, 1)
  T = knockoff.threshold(W, q, offset)   # threshold

  selected = (1:p)[W >= T]

}







