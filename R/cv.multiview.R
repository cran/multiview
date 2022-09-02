#' Perform k-fold cross-validation for cooperative learning
#'
#' Does k-fold cross-validation (CV) for multiview and produces a CV curve.
#' 
#' The current code can be slow for "large" data sets, e.g. when the
#' number of features is larger than 1000.  It can be helpful to see
#' the progress of multiview as it runs; to do this, set trace.it = 1
#' in the call to multiview or cv.multiview.  With this, multiview
#' prints out its progress along the way.  One can also pre-filter the
#' features to a smaller set, using the exclude option, with a filter
#' function.
#' 
#' If there are missing values in the feature matrices: 
#' we recommend that you center the columns of each feature matrix, and then fill in the missing values with 0.
#' 
#' For example, \cr
#' `x <- scale(x,TRUE,FALSE)` \cr
#' `x[is.na(x)] <- 0` \cr
#' `z <- scale(z,TRUE,FALSE)` \cr
#' `z[is.na(z)] <- 0` \cr
#' 
#' Then run multiview in the usual way. It will exploit the assumed shared latent factors
#' to make efficient use of the available data.
#'
#' The function runs `multiview` `nfolds+1` times; the first to get the
#' `lambda` sequence, and then the remainder to compute the fit with each
#' of the folds omitted. The error is accumulated, and the average error and
#' standard deviation over the folds is computed.  Note that `cv.multiview`
#' does NOT search for values for `rho`. A specific value should be
#' supplied, else `rho=0` is assumed by default. If users would like to
#' cross-validate `rho` as well, they should call `cv.multiview` with
#' a pre-computed vector `foldid`, and then use this same fold vector in
#' separate calls to `cv.multiview` with different values of `rho`.
#'
#' @inheritParams multiview
#' @param weights Observation weights; defaults to 1 per observation
#' @param offset Offset vector (matrix) as in `multiview`
#' @param type.measure loss to use for cross-validation. Currently
#'   five options, not all available for all models. The default is
#'   `type.measure="deviance"`, which uses squared-error for gaussian
#'   models (a.k.a `type.measure="mse"` there), deviance for logistic
#'   and poisson regression, and partial-likelihood for the Cox model.
#'   `type.measure="class"` applies to binomial and multinomial
#'   logistic regression only, and gives misclassification error.
#'   `type.measure="auc"` is for two-class logistic regression only,
#'   and gives area under the ROC curve. `type.measure="mse"` or
#'   `type.measure="mae"` (mean absolute error) can be used by all
#'   models except the `"cox"`; they measure the deviation from the
#'   fitted mean to the response.  `type.measure="C"` is Harrel's
#'   concordance measure, only available for `cox` models.
#' @param nfolds number of folds - default is 10. Although `nfolds`
#'   can be as large as the sample size (leave-one-out CV), it is not
#'   recommended for large datasets. Smallest value allowable is
#'   `nfolds=3`
#' @param foldid an optional vector of values between 1 and `nfold`
#'   identifying what fold each observation is in. If supplied,
#'   `nfold` can be missing.
#' @param mvlambda A user supplied `lambda` sequence, default
#'   `NULL`. Typical usage is to have the program compute its own
#'   `mvlambda` sequence. This sequence, in general, is different from
#'   that used in the [glmnet::glmnet()] call (named `lambda`). Note
#'   that this is done for the full model (master sequence), and
#'   separately for each fold.  The fits are then aligned using the
#'   glmnet lambda sequence associated with the master sequence (see
#'   the `alignment` argument for additional details). Adapting
#'   `mvlambda` for each fold leads to better convergence. When
#'   `mvlambda` is supplied, the same sequence is used everywhere, but
#'   in some GLMs can lead to convergence issues.
#' @param alignment This is an experimental argument, designed to fix
#'   the problems users were having with CV, with possible values
#'   `"lambda"` (the default) else `"fraction"`. With `"lambda"` the
#'   `lambda` values from the master fit (on all the data) are used to
#'   line up the predictions from each of the folds. In some cases
#'   this can give strange values, since the effective `lambda` values
#'   in each fold could be quite different. With `"fraction"` we line
#'   up the predictions in each fold according to the fraction of
#'   progress along the regularization. If in the call a `lambda`
#'   argument is also provided, `alignment="fraction"` is ignored
#'   (with a warning).
#' @param grouped This is an experimental argument, with default
#'   `TRUE`, and can be ignored by most users. For all models except
#'   the `"cox"`, this refers to computing `nfolds` separate
#'   statistics, and then using their mean and estimated standard
#'   error to describe the CV curve. If `grouped=FALSE`, an error
#'   matrix is built up at the observation level from the predictions
#'   from the `nfold` fits, and then summarized (does not apply to
#'   `type.measure="auc"`). For the `"cox"` family, `grouped=TRUE`
#'   obtains the CV partial likelihood for the Kth fold by
#'   \emph{subtraction}; by subtracting the log partial likelihood
#'   evaluated on the full dataset from that evaluated on the on the
#'   (K-1)/K dataset. This makes more efficient use of risk sets. With
#'   `grouped=FALSE` the log partial likelihood is computed only on
#'   the Kth fold
#' @param keep If `keep=TRUE`, a \emph{prevalidated} array is returned
#'   containing fitted values for each observation and each value of
#'   `lambda`. This means these fits are computed with this
#'   observation and the rest of its fold omitted. The `foldid` vector
#'   is also returned.  Default is keep=FALSE.
#' @param trace.it If `trace.it=1`, then progress bars are displayed;
#'   useful for big models that take a long time to fit.
#' @param \dots Other arguments that can be passed to `multiview`
#' @return an object of class `"cv.multiview"` is returned, which is a
#'   list with the ingredients of the cross-validation
#'   fit. \item{lambda}{the values of `lambda` used in the fits.}
#'   \item{cvm}{The mean cross-validated error - a vector of length
#'   `length(lambda)`.} \item{cvsd}{estimate of standard error of
#'   `cvm`.} \item{cvup}{upper curve = `cvm+cvsd`.} \item{cvlo}{lower
#'   curve = `cvm-cvsd`.} \item{nzero}{number of non-zero coefficients
#'   at each `lambda`.} \item{name}{a text string indicating type of
#'   measure (for plotting purposes).} \item{multiview.fit}{a fitted
#'   multiview object for the full data.} \item{lambda.min}{value of
#'   `lambda` that gives minimum `cvm`.} \item{lambda.1se}{largest
#'   value of `lambda` such that error is within 1 standard error of
#'   the minimum.} \item{fit.preval}{if `keep=TRUE`, this is the array
#'   of prevalidated fits. Some entries can be `NA`, if that and
#'   subsequent values of `lambda` are not reached for that fold}
#'   \item{foldid}{if `keep=TRUE`, the fold assignments used}
#'   \item{index}{a one column matrix with the indices of `lambda.min`
#'   and `lambda.1se` in the sequence of coefficients, fits etc.}
#'
#' @examples
#' # Gaussian
#' # Generate data based on a factor model
#' set.seed(1)
#' x = matrix(rnorm(100*20), 100, 20)
#' z = matrix(rnorm(100*20), 100, 20)
#' U = matrix(rnorm(100*5), 100, 5)
#' for (m in seq(5)){
#'     u = rnorm(100)
#'     x[, m] = x[, m] + u
#'     z[, m] = z[, m] + u
#'     U[, m] = U[, m] + u}
#' x = scale(x, center = TRUE, scale = FALSE)
#' z = scale(z, center = TRUE, scale = FALSE)
#' beta_U = c(rep(0.1, 5))
#' y = U %*% beta_U + 0.1 * rnorm(100)
#' fit1 = cv.multiview(list(x=x,z=z), y, rho = 0.3)
#' 
#' # plot the cross-validation curve
#' plot(fit1)
#' 
#' # extract coefficients
#' coef(fit1, s="lambda.min")
#' 
#' # extract ordered coefficients
#' coef_ordered(fit1, s="lambda.min")
#' 
#' # make predictions
#' predict(fit1, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min")
#' 
#' # Binomial
#' \donttest{
#' by = 1 * (y > median(y)) 
#' fit2 = cv.multiview(list(x=x,z=z), by, family = binomial(), rho = 0.9)
#' predict(fit2, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min", type = "response")
#' plot(fit2)
#' coef(fit2, s="lambda.min")
#' coef_ordered(fit2, s="lambda.min")
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = cv.multiview(list(x=x,z=z), py, family = poisson(), rho = 0.6)
#' predict(fit3, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min", type = "response") 
#' plot(fit3)
#' coef(fit3, s="lambda.min")
#' coef_ordered(fit3, s="lambda.min")
#' }
#' @importFrom glmnet glmnet
#' @importFrom stats median
#' @export
cv.multiview <- function(x_list, y, family = gaussian(), rho = 0, weights = NULL, offset = NULL, mvlambda = NULL,
                         type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
                         nfolds = 10, foldid = NULL,  alignment = c("lambda", "fraction"),
                         grouped = TRUE, keep = FALSE, trace.it = 0, ...) {

  type.measure <-  match.arg(type.measure)
  alignment <- match.arg(alignment)
  if (!is.null(mvlambda) && length(mvlambda) < 2)
    stop("Need more than one value of lambda for cv.multiview")
  if (!is.null(mvlambda) && alignment == "fraction"){
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment <- "lambda"
  }
  m  <- length(x_list)
  N <- nrow(x_list[[1L]])

  y <- drop(y)

  cv.call <- multiview.call <- match.call(expand.dots = TRUE)
  which <- match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(multiview.call), FALSE)
  if (any(which)) {
    multiview.call <- multiview.call[-which]
  }
  multiview.call[[1]] <- as.name("multiview")

  if (multiview.control()$itrace) {
    trace.it <- 1
  } else {
    if(trace.it) {
      multiview.control(itrace = 1)
      on.exit(multiview.control(itrace = 0))
    }
  }

  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else {
    nfolds <- max(foldid)
  }
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }

  p <- do.call(sum, lapply(x_list, ncol))

  if (is.null(weights)) {
    weights  <- rep(1.0, N)
  }

  if (trace.it) cat("Training\n")
  multiview.object <- multiview(x_list = x_list, y = y, family = family, rho = rho,
                                weights = weights, offset = offset,
                                mvlambda = mvlambda,
                                trace.it = trace.it, ...)
  multiview.object$call <- multiview.call
  lognet_class <- multiview.object$family$family == "binomial"
  if (lognet_class) {
    subclass <- "lognet"
  } else {
    subclass <- class(multiview.object)[[1L]]
  }
  type.measure <- cvtype(type.measure, subclass)
  is.offset <- multiview.object$offset
  ## Next line is commented out so each call generates its own lambda sequence
  ## lambda=glmnet.object$lambda
  #dlist <- list()
  if (inherits(multiview.object, "multnet") && !multiview.object$grouped) {
    nz <- predict(multiview.object, type = "nonzero")
    nz <- sapply(nz, function(x) sapply(x, length))
    nz <- ceiling(apply(nz, 1, median))
  } else {
    nz <-  sapply(predict(multiview.object, type = "nonzero"),
                  length)
  }
  outlist <- as.list(seq(nfolds))

  for (i in seq(nfolds)){
    if (trace.it) cat(sprintf("Fold: %d/%d\n", i, nfolds))
    which <-  (foldid == i)
    x_sub_list <- lapply(x_list, function(x) x[!which, , drop = FALSE])
    y_sub <- y[!which, drop = FALSE]
    if (is.null(weights)) {
      weights_sub <- NULL
    } else {
      weights_sub <- weights[!which]
    }
    if (is.null(offset)) {
      offset_sub  <- NULL
    } else {
      offset_sub <- offset[!which, drop = FALSE]
    }
    #dlist[[i]]  <- list(x = x_sub_list, y = y_sub)
    outlist[[i]] <- multiview(x_list = x_sub_list, y = y_sub, family = family, rho = rho,
                              weights = weights_sub, offset = offset_sub,
                              mvlambda = mvlambda, trace.it = trace.it, ...)
  }
  #saveRDS(dlist, "m_dlist.RDS")
  #saveRDS(outlist, "m_outlist.RDS")
  lambda <- multiview.object$lambda
  class(outlist) <- paste0(subclass, "list")

  predmat <- build_predmat(outlist, lambda, x_list, offset, foldid, alignment, y = y, weights = weights,
                           grouped = grouped, type.measure = type.measure, family = family(multiview.object)
                           )
  ### we include type.measure for the special case of coxnet with the deviance vs C-index discrepancy
  ### family is included for the new GLM crowd
  ### Next we compute the measures
  if(subclass == "multiview") attr(predmat, "family") <- multiview.object$family
  #fun <- paste("cv", subclass, sep = ".")
  #cvstuff <- do.call(fun, list(predmat, y, type.measure, weights, foldid, grouped))
  ## The above (from glmnet) needs to be modified for our multiview case. We merely call
  ## cv.lognet in case of binomial.
  if (lognet_class) {
    cvstuff <- cv.lognet(predmat, y, type.measure, weights, foldid, grouped)
  } else {
    #cvstuff <- do.call(cv.glmnetfit, list(predmat, y, type.measure, weights, foldid, grouped))
    cvstuff <- cv.glmnetfit(predmat, y, type.measure, weights, foldid, grouped)
  }

  grouped <- cvstuff$grouped
  if ((N / nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.multiview, since < 3 observations per fold",
            call. = FALSE)
    grouped <- FALSE
  }

  out <- cvstats(cvstuff, foldid, nfolds, lambda, nz, grouped)
  cvname <- names(cvstuff$type.measure)
  names(cvname) <- cvstuff$type.measure# to be compatible with earlier version; silly, I know
  out <- c(out, list(call = cv.call, name = cvname, multiview.fit = multiview.object))
  if (keep) {
    out <- c(out, list(fit.preval = predmat, foldid = foldid))
  }
  lamin <- with(out, getOptcv.glmnet(lambda, cvm, cvsd, cvname))
  obj <- c(out, as.list(lamin))
  class(obj) <- c("cv.multiview", "cv.glmnet")
  obj
}

#' Extract coefficients from a cv.multiview object
#'
#' @param object Fitted `"cv.multiview"` object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the value `s="lambda.1se"`
#'   stored on the CV `object`. Alternatively `s="lambda.min"` can be
#'   used. If `s` is numeric, it is taken as the value(s) of `lambda`
#'   to be used. (For historical reasons we use the symbol 's' rather
#'   than 'lambda' to reference this parameter.)
#' @return the matrix of coefficients for specified lambda.
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100*20), 100, 20)
#' z = matrix(rnorm(100*20), 100, 20)
#' U = matrix(rnorm(100*5), 100, 5)
#' for (m in seq(5)){
#'     u = rnorm(100)
#'     x[, m] = x[, m] + u
#'     z[, m] = z[, m] + u
#'     U[, m] = U[, m] + u}
#' x = scale(x, center = TRUE, scale = FALSE)
#' z = scale(z, center = TRUE, scale = FALSE)
#' beta_U = c(rep(0.1, 5))
#' y = U %*% beta_U + 0.1 * rnorm(100)
#' fit1 = cv.multiview(list(x=x,z=z), y, rho = 0.3)
#' coef(fit1, s="lambda.min")
#' 
#' # Binomial
#' \donttest{
#' by = 1 * (y > median(y)) 
#' fit2 = cv.multiview(list(x=x,z=z), by, family = binomial(), rho = 0.9)
#' coef(fit2, s="lambda.min")
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = cv.multiview(list(x=x,z=z), py, family = poisson(), rho = 0.6)
#' coef(fit3, s="lambda.min")
#' }
#' @importFrom stats coef 
#' @method coef cv.multiview
#' @inheritParams predict.multiview
#' @export
coef.cv.multiview <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  if(is.numeric(s))
    lambda <- s
  else
    if(is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    }
  else
    stop("Invalid form for s")
  coef(object$multiview.fit, s = lambda, ...)
}

#' Extract an ordered list of standardized coefficients from a cv.multiview object
#'
#' This function extracts a ranked list of coefficients after the coefficients are
#' standardized by the standard deviation of the corresponding features. The ranking
#' is based on the magnitude of the standardized coefficients. It also outputs
#' the data view to which each coefficient belongs.
#'
#' The output table shows from left to right the data view each coefficient comes from,
#' the column index of the feature in the corresponding data view, the coefficient
#' after being standardized by the standard deviation of the corresponding feature,
#' and the original fitted coefficient.
#'
#' @param object Fitted `"cv.multiview"` object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#' predictions are required. Default is the value `s="lambda.1se"` stored
#' on the CV `object`. Alternatively `s="lambda.min"` can be used. If
#' `s` is numeric, it is taken as the value(s) of `lambda` to be
#' used. (For historical reasons we use the symbol 's' rather than 'lambda' to
#' reference this parameter.)
#' @return data frame of consisting of view name, view column,
#'   coefficient and standardized coefficient ordered by rank of
#'   standardized coefficient.
#'
#' @examples
#' set.seed(1)
#' x = matrix(rnorm(100*20), 100, 20)
#' z = matrix(rnorm(100*20), 100, 20)
#' U = matrix(rnorm(100*5), 100, 5)
#' for (m in seq(5)){
#'     u = rnorm(100)
#'     x[, m] = x[, m] + u
#'     z[, m] = z[, m] + u
#'     U[, m] = U[, m] + u}
#' x = scale(x, center = TRUE, scale = FALSE)
#' z = scale(z, center = TRUE, scale = FALSE)
#' beta_U = c(rep(0.1, 5))
#' y = U %*% beta_U + 0.1 * rnorm(100)
#' fit1 = cv.multiview(list(x=x,z=z), y, rho = 0.3)
#' coef_ordered(fit1, s="lambda.min")
#' 
#' # Binomial
#' \donttest{
#' by = 1 * (y > median(y)) 
#' fit2 = cv.multiview(list(x=x,z=z), by, family = binomial(), rho = 0.9)
#' coef_ordered(fit2, s="lambda.min")
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = cv.multiview(list(x=x,z=z), py, family = poisson(), rho = 0.6)
#' coef_ordered(fit3, s="lambda.min")
#' }
#' @method coef_ordered cv.multiview
#' @inheritParams predict.multiview
#' @export
coef_ordered.cv.multiview <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  if (length(s) != 1) {
    stop("s has to be specified as a single value")
  }
  if(is.numeric(s))
    lambda <- s
  else
    if(is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    }
  else
    stop("Invalid form for s")
  coef_ordered(object$multiview.fit, s = lambda, ...)
}




