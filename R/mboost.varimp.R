
stat.mboost_varimp2 <- function(Xaug, y, max.mstop = 100, bl = c("bbs", "bols", "btree"),
                                cv.fold = 5, family = Gaussian()){
  library(mboost)
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
    Xaug$y = y
  }
  model = mboost(y ~ ., data = Xaug, control = boost_control(mstop = max.mstop), baselearner = bl, family = family)
  cv10f = cv(model.weights(model), type="kfold", B = cv.fold)
  cvm = cvrisk(model, folds = cv10f, papply = mclapply)
  best.model = model[mstop(cvm)]
  Z = as.numeric(varimp(best.model))
}
