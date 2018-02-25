lasso_coef_diff <- function (Xaug, y, family = Gaussian())
{
  check.dist <- function(x) grepl(x,family@name)
  if (any(unlist(lapply(c("Ada","Binomial","AUC"),check.dist)))) response.type = "binary"
  else if (any(unlist(lapply(c("Squared Error","Huber","Absolute Err"),check.dist)))) response.type = "continuous"
  else if (any(unlist(lapply(c("Poisson","Gamma","Multinomial"),check.dist)))) response.type = "discrete"
  else if (any(unlist(lapply(c("Cox","Weibull","Log Logistic","Lognormal","Gehan","Concordance Probability"),check.dist))))  response.type = "survival"
  else stop("unknown family")

  if (response.type %in% c("continuous","discrete")){
    y = as.vector(y)
    fam = "gaussian"
  }
  else if (response.type == "binary"){
    y = as.factor(y)
    fam = "binomial"
  }
  else fam = "cox"

  if (!requireNamespace("glmnet", quietly = T))
    stop("glmnet is not installed", call. = F)
  parallel = T
  if (!requireNamespace("doMC", quietly = T)) {
    warning("doMC is not installed. Without parallelization, the statistics will be slower to compute",
            call. = F, immediate. = T)
    parallel = F
  }
  if (!requireNamespace("parallel", quietly = T)) {
    warning("parallel is not installed. Without parallelization, the statistics will be slower to compute.",
            call. = F, immediate. = T)
    parallel = F
  }
  if (parallel) {
    ncores = parallel::detectCores(all.tests = TRUE, logical = TRUE)
    cores = min(2, ncores)
    if (cores > 1) {
      doMC::registerDoMC(cores = cores)
      parallel = TRUE
    }
    else {
      parallel = FALSE
    }
  }
  Z = cv_coeffs_glmnet(Xaug, y, family = fam,
                       parallel = parallel)
  Z = abs(Z)
}


#' @keywords internal
cv_coeffs_glmnet <- function(X, y, family = "gaussian", parallel=T) {
  # Standardize variables
  X = normc(X)

  n = nrow(X); p = ncol(X)
  nlambda = 100

    if( identical(family, "gaussian") ) {
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }

  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda=lambda, family = family,
                                     standardize=F,standardize.response=F, parallel=parallel)

  coef(cv.glmnet.fit, s = "lambda.min")[2:(p+1)]
}

normc <- function (X)
{
  X.scaled = scale(X, center = TRUE, scale = FALSE)
  X.scaled = scale(X.scaled, center = FALSE, scale = sqrt(colSums(X.scaled^2)))
  X.scaled[, ]
}

