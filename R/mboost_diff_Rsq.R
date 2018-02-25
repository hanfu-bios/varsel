
stat.mboost_diff_Rsq <- function(Xaug, y, max.mstop = 100, bl = c("bbs", "bols", "btree"),
                                cv.fold = 5, family = Gaussian()){

  p = ncol(Xaug) / 2
  Xaug = as.data.frame(Xaug)
  check.dist <- function(x) grepl(x,family@name)
  if (any(unlist(lapply(c("Ada","Binomial","AUC"),check.dist)))) response.type = "binary"
  else if (any(unlist(lapply(c("Squared Error","Huber","Absolute Err"),check.dist)))) response.type = "continuous"
  else if (any(unlist(lapply(c("Poisson","Gamma","Multinomial"),check.dist)))) response.type = "discrete"
  else if (any(unlist(lapply(c("Cox","Weibull","Log Logistic","Lognormal","Gehan","Concordance Probability"),check.dist))))  response.type = "survival"
  else stop("unknown family")

  if (response.type %in% c("continuous","discrete")){
    Xaug$y = as.vector(y)
  }
  else if (response.type == "binary"){
    Xaug$y = as.factor(y)
  }
  else {   # survival
    error("R-square for survival data still under development")
  }
  model = mboost(y ~ ., data = Xaug, control = boost_control(mstop = max.mstop), baselearner = bl, family = family)
  cv10f = cv(model.weights(model), type="kfold", B = cv.fold)
  cvm = cvrisk(model, folds = cv10f, papply = mclapply)

  best.mstop = mstop(cvm)
  best.fit = cv.mboost(Xaug, y, best.mstop, family = family, baselearner = bl, cv.fold = cv.fold)
  mboost.mse = mean((y-best.fit)^2)
  Rsq = mse2Rsq(mboost.mse, y)
  mse_j = NULL
  Rsq_j = NULL
  for (j in 1:(2*p)){
    best.pred1 = cv.mboost(Xaug[,-j], y, best.mstop, family = family, baselearner = bl, cv.fold = cv.fold)
    mse_j[j] = mean((best.pred1 - y)^2)
    Rsq_j[j] = mse2Rsq(mse_j[j], y)
  }
  Z = abs(Rsq - Rsq_j)

}

mse2Rsq <- function(mse, y) 1 - mse*n/sum((y-mean(y))^2)

cv.mboost <- function(Xaug, y, mstop, family, baselearner = "bbs", cv.fold = 5){
  folds_i = sample(rep(1:cv.fold, length.out = n))
  X.pred2 = rep(0,n)
  for (k2 in 1:cv.fold) {
    test_i = which(folds_i == k2)
    Xtr = Xaug[-test_i, ]
    Xt = Xaug[test_i, ]
    fit1 = mboost(y ~ ., data = Xtr, control = boost_control(mstop = mstop), baselearner = baselearner, family = family)
    X.pred2[test_i] = predict(fit1, newdata = subset(Xt,select=-y), type = "response")
  }
  X.pred2
}
