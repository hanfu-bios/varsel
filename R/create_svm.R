create.svm <- function(X){
  Xaug = X
  n = nrow(X)
  p = ncol(X)
  X_k = matrix(0, n, p)
  ncost = 20
  cost = 2^(seq(-5, 12, length.out=ncost))
  # cost = seq(1, 400, length.out=ncost)
  gamma = 1
  # gamma = 0.0001
  # epsilon = 0.1
  # nfolds = 10
  # tolerance = 0.003
  for (j in 1:p){
    X.pred = matrix(0, n, ncost)
    error = NULL
    for (c in 1:ncost){
      # X.model = svm(Xaug[,j] ~ ., data = Xaug[,-j], cost = cost[c], cross = nfolds, gamma = gamma,
      #               epsilon = epsilon, scale = F)
      # error[c] = X.model$tot.MSE
      # X.pred[,c] = predict(X.model)
      X.pred[,c] = cv.svm(Xaug[,-j], Xaug[,j], cost[c], gamma)
      error[c] = mean((X.pred[,c]-Xaug[,j])^2)
    }
    # best = which(error < mean(error[(ncost-20):ncost])+tolerance)[1]
    best = which.min(error)
    X.best.pred = X.pred[,best]
    err = Xaug[,j] - X.best.pred
    error_permute = err[sample(n)]
    X_k[,j] = X.best.pred + error_permute
    Xaug = cbind(Xaug, X_k[,j])
  }
  return(X_k)
}


cv.svm <- function(x, y, cost, gamma){
  # set.seed(1234)
  n_folds = 5
  folds_i = sample(rep(1:n_folds, length.out = n))
  X.pred2 = rep(0,n)
  # train.error = NULL
  # test.error = NULL
  # error.all = NULL
  for (k2 in 1:n_folds) {
    test_i = which(folds_i == k2)
    Xtr = x[-test_i, ]
    Xt = x[test_i, ]
    ytr = y[-test_i]
    yt = y[test_i]
    X.model2 = svm(x = Xtr, y = ytr, cost = cost, gamma = gamma)
    # X.pred1 = predict(X.model2,Xtr)
    # train.error[k2] = mean((Xaugtr[,1]-X.pred1)^2)
    X.pred2[test_i] = predict(X.model2, Xt)
    # error.all[k2] = mean((X.pred2[test_i]-yt)^2)
  }
  # list = list(pred=X.pred2, error.all = error.all)
  return(X.pred2)
  # list
}
