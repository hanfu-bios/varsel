
stat.mboost_varimp2 <- function(Xaug, y, response.type, max.mstop, bl, cv.fold, family){
  library(mboost)
  Xaug = as.data.frame(Xaug)
  if (response.type == "continuous"){
    Xaug$y = as.vector(y)
    if(is.null(family)) family = Gaussian()
  }
  else if (response.type == "binomial"){
    Xaug$y = as.factor(y)
    if(is.null(family)) family = Binomial()
  }
  else {   # survival
    Xaug$y = y
    if(is.null(family)) family = CoxPH()
  }
  model = mboost(y ~ ., data = Xaug, control = boost_control(mstop = max.mstop), baselearner = bl, family = family)
  cv10f = cv(model.weights(model), type="kfold", B = cv.fold)
  cvm = cvrisk(model, folds = cv10f, papply = mclapply)
  best.model = model[mstop(cvm)]
  Z = as.numeric(varimp(best.model))
}
