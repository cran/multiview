#' Fit a GLM with elastic net regularization for a path of lambda values
#'
#' Fit a generalized linear model via penalized maximum likelihood for a path of
#' lambda values. Can deal with any GLM family.
#'
#' `multiview.path` solves the elastic net problem for a path of lambda values.
#' It generalizes `multiview::multiview` in that it works for any GLM family.
#'
#' Sometimes the sequence is truncated before `nlam` values of lambda
#' have been used. This happens when `multiview.path` detects that the decrease
#' in deviance is marginal (i.e. we are near a saturated fit).
#'
#' @inheritParams multiview
#' @param x the `cbind`ed matrices in `x_list`
#' @param nvars the number of variables (total)
#' @param nobs the number of observations
#' @param xm the column means vector (could be zeros if `standardize = FALSE`)
#' @param xs the column std dev vector (could be 1s if `standardize = FALSE`)
#' @param control the multiview control object
#' @param vp the variable penalities (processed)
#' @param vnames the variable names
#' @param is.offset a flag indicating if offset is supplied or not
#' @param user_lambda a flag indicating if user supplied the lambda sequence
#' @param start_val the result of first call to `get_start`
#' @return An object with class `"multiview"` `"glmnetfit"` and `"glmnet"`
#' \item{a0}{Intercept sequence of length `length(lambda)`.}
#' \item{beta}{A `nvars x length(lambda)` matrix of coefficients, stored in
#' sparse matrix format.}
#' \item{df}{The number of nonzero coefficients for each value of lambda.}
#' \item{dim}{Dimension of coefficient matrix.}
#' \item{lambda}{The actual sequence of lambda values used. When alpha=0, the
#' largest lambda reported does not quite give the zero coefficients reported
#' (lambda=inf would in principle). Instead, the largest lambda for alpha=0.001
#' is used, and the sequence of lambda values is derived from this.}
#' \item{lambda}{The sequence of lambda values}
#' \item{mvlambda}{The corresponding sequence of multiview lambda values}
#' \item{dev.ratio}{The fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
#' for the saturated model (a model with a free parameter per observation).
#' Hence dev.ratio=1-dev/nulldev.}
#' \item{nulldev}{Null deviance (per observation). This is defined to be
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data summed over all lambda values.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{A logical variable indicating whether an offset was included
#' in the model.}
#' \item{call}{The call that produced this object.}
#' \item{family}{Family used for the model.}
#' \item{nobs}{Number of observations.}
#'
#' @import methods
#' @importFrom Matrix Matrix
# @examples
# set.seed(1)
# x <- matrix(rnorm(100 * 20), nrow = 100)
# z <- matrix(rnorm(100 * 20), nrow = 100)
# y <- ifelse(rnorm(100) > 0, 1, 0)
#
# # binomial with probit link
# fit1 <- multiview:::multiview.path(list(x, z), y, family = binomial(link = "probit"),
#                                    x = cbind(x, z))
multiview.path <- function(x_list, y, rho = 0, weights = NULL,
                           lambda, nlambda, user_lambda = FALSE,
                           alpha = 1.0, offset = NULL, family = gaussian(),
                           standardize = TRUE, intercept = TRUE, thresh = 1e-7, maxit = 100000,
                           penalty.factor = rep(1.0, nvars), exclude = integer(0), lower.limits = -Inf,
                           upper.limits = Inf, trace.it = 0, x, nvars, nobs, xm, xs, control, vp, vnames, start_val, is.offset) {

## multiview.path <- function(x_list, y, rho = 0, weights = NULL, mvlambda = NULL, nlambda = 100,
##                            lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001),
##                            alpha = 1.0, offset = NULL, family = gaussian(),
##                            standardize = TRUE, intercept = TRUE, thresh = 1e-7, maxit = 100000,
##                            penalty.factor = rep(1.0, nvars), exclude = integer(0), lower.limits = -Inf,
##                            upper.limits = Inf, trace.it = 0, x) {

  ##   ### Check on family argument
  ## if (is.function(family)) {
  ##   family <- family()
  ## }
  
  ## ## Make lambda the multiview lambda
  ## ## Allows use of mvlambda parameter for idempotence.
  ## lambda <- mvlambda
  
  ## this.call <- match.call()

  ## ## ## Prepare to reuse glmnetFlex code
  ## ## x <- do.call(cbind, x_list)

  ## ## ## We need the std devs for other purposes, so we compute it
  ## ## xsd <- apply(x, 2, sd)

  ## ### Prepare all the generic arguments
  ## ## if (alpha > 1) {
  ## ##   warning("alpha > 1; set to 1")
  ## ##   alpha = 1
  ## ## } else if (alpha < 0) {
  ## ##   warning("alpha < 0; set to 0")
  ## ##   alpha = 0
  ## ## }
  ## alpha = as.double(alpha)

  ## np = dim(x)
  ## #if(is.null(np) || (np[2] <= 1)) stop("x should be a matrix with 2 or more columns")
  ## nobs = as.integer(np[1]); nvars = as.integer(np[2])

  ## # get feature variable names
  ## vnames=colnames(x)
  ## #if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")

  ## # check weights
  ## if(is.null(weights)) weights = rep(1,nobs)
  ## else if (length(weights) != nobs)
  ##   stop(paste("Number of elements in weights (",length(weights),
  ##              ") not equal to the number of rows of x (",nobs,")",sep=""))
  ## weights <- as.double(weights)

  ## ## initialize from family function. Makes y a vector in case of binomial, and possibly changes weights
  ## ## Expects nobs to be defined, and creates n and mustart (neither used here)
  ## ## Some cases expect to see things, so we set it up just to make it work
  ## etastart=0;mustart=NULL;start=NULL
  ## eval(family$initialize)
  ## ##
  ## ## Just in case this was not done in initialize()
  ## y <- drop(y)  # we don't like matrix responses

  ## is.offset <- !(is.null(offset))
  ## if (is.offset == FALSE) {
  ##   offset <- as.double(y * 0) #keeps the shape of y
  ## }
  ## # infinite penalty factor vars are excluded
  ## if(any(penalty.factor == Inf)) {
  ##   exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
  ##   exclude = sort(unique(exclude))
  ## }

  ## ## Compute weighted mean and variance of columns of x, sensitive to sparse matrix
  ## ## needed to detect constant columns below, and later if standarization
  ## meansd <- weighted_mean_sd(x, weights)

  ## ## look for constant variables, and if any, then add to exclude
  ## const_vars <- meansd$sd == 0
  ## nzvar <- setdiff(which(!const_vars), exclude)
  ## # if all the non-excluded variables have zero variance, throw error
  ## if (length(nzvar) == 0) stop("All used predictors have zero variance")

  ## ## if any constant vars, add to exclude
  ## if(any(const_vars)) {
  ##   exclude <- sort(unique(c(which(const_vars),exclude)))
  ##   meansd$sd[const_vars] <- 1.0 ## we divide later, and do not want bad numbers
  ## }
  ## if(length(exclude) > 0) {
  ##   jd = match(exclude, seq(nvars), 0)
  ##   if(!all(jd > 0)) stop ("Some excluded variables out of range")
  ##   penalty.factor[jd] = 1 # ow can change lambda sequence
  ## }
  ## # check and standardize penalty factors (to sum to nvars)
  ## vp = pmax(0, penalty.factor)
  ## if (max(vp) <= 0) stop("All penalty factors are <= 0")
  ## vp = as.double(vp * nvars / sum(vp))


  ## ### check on limits
  ## control <- multiview.control()
  ## if (thresh >= control$epsnr)
  ##   warning("thresh should be smaller than multiview.control()$epsnr",
  ##           call. = FALSE)

  ## if(any(lower.limits > 0)){ stop("Lower limits should be non-positive") }
  ## if(any(upper.limits < 0)){ stop("Upper limits should be non-negative") }
  ## lower.limits[lower.limits == -Inf] = -control$big
  ## upper.limits[upper.limits == Inf] = control$big
  ## if (length(lower.limits) < nvars) {
  ##   if(length(lower.limits) == 1) lower.limits = rep(lower.limits, nvars) else
  ##                                                                           stop("Require length 1 or nvars lower.limits")
  ## } else lower.limits = lower.limits[seq(nvars)]
  ## if (length(upper.limits) < nvars) {
  ##   if(length(upper.limits) == 1) upper.limits = rep(upper.limits, nvars) else
  ##                                                                           stop("Require length 1 or nvars upper.limits")
  ## } else upper.limits = upper.limits[seq(nvars)]

  ## if (any(lower.limits == 0) || any(upper.limits == 0)) {
  ##   ###Bounds of zero can mess with the lambda sequence and fdev;
  ##   ###ie nothing happens and if fdev is not zero, the path can stop
  ##   fdev <- multiview.control()$fdev
  ##   if(fdev!= 0) {
  ##     multiview.control(fdev = 0)
  ##     on.exit(multiview.control(fdev = fdev))
  ##   }
  ## }
  ## ### end check on limits
  ## ### end preparation of generic arguments

  ## # standardize x if necessary
  ## ## if (intercept) {
  ## ##   xm <- meansd$mean
  ## ## } else {
  ## ##   xm <- rep(0.0, times = nvars)
  ## ## }
  ## ## We handle intercept ourselves!
  ## xm <- rep(0.0, times = nvars)
  ## if (standardize) {
  ##   xs <- meansd$sd
  ## } else {
  ##   xs <- rep(1.0, times = nvars)
  ## }
  ## if (!inherits(x, "sparseMatrix")) {
  ##   x <- scale(x, xm, xs)
  ## } else {
  ##   attr(x, "xm") <- xm
  ##   attr(x, "xs") <- xs
  ## }
  ## lower.limits <- lower.limits * xs
  ## upper.limits <- upper.limits * xs

  ## # get null deviance and lambda max
  ## start_val <- get_start(x, y, weights, family, intercept, is.offset,
  ##                        offset, exclude, vp, alpha)

  ## # work out lambda values
  ## nlam = as.integer(nlambda)
  ## user_lambda = FALSE   # did user provide their own lambda values?
  ## if (is.null(lambda)) {
  ##   if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")

  ##   # compute lambda max: to add code here
  ##   lambda_max <- start_val$lambda_max

  ##   # compute lambda sequence
  ##   ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
  ##                   length.out = nlam))
  ## } else {  # user provided lambda values
  ##   user_lambda = TRUE
  ##   if (any(lambda < 0)) stop("lambdas should be non-negative")
  ##   ulam = as.double(rev(sort(lambda)))
  ##   nlam = as.integer(length(lambda))
  ## }


  ### NOTA BENE
  ## Up to this everything is already set up now because of standardization.
  ## END OF NB

  # start progress bar
  if (trace.it == 1) pb <- utils::txtProgressBar(min = 0, max = nlambda, style = 3)

  
  glambda <- rep(1.0, nlambda) # the actual glmnet lambda sequence, initially the scale factor
  a0 <- rep(NA, length = nlambda)
  beta <- matrix(0, nrow = nvars, ncol = nlambda)
  dev.ratio <- rep(NA, length = nlambda)
  fit <- NULL
  mnl <- min(nlambda, control$mnlam)
  cur_lambda <- lambda
  cur_lambda[1] <- if(user_lambda) lambda[1] else control$big
  for (k in 1:nlambda) {
    # get the correct lambda value to fit
    ## if (k > 1) {
    ##   cur_lambda <- ulam[k]
    ## } else {
    ##   cur_lambda <- if(user_lambda) ulam[k] else control$big
    ## }
    ## effective_lambda <- cur_lambda[k]
    if (trace.it == 2) cat("Fitting lambda index", k, ":", lambda[k], fill = TRUE)
    fit <- multiview.fit(x_list = x_list,
                         x = x,
                         y = y,
                         rho = rho,
                         #weights = weights / sum(weights),
                         weights = weights,
                         lambda = cur_lambda[k],
                         alpha = alpha,
                         offset = offset,
                         family = family,
                         intercept = intercept,
                         thresh = thresh,
                         maxit = maxit,
                         penalty.factor = vp,
                         exclude = exclude,
                         lower.limits = lower.limits,
                         upper.limits = upper.limits,
                         warm = fit,
                         from.multiview.path = TRUE,
                         save.fit = TRUE,
                         trace.it = trace.it,
                         user_lambda = user_lambda)
    if (trace.it == 1) utils::setTxtProgressBar(pb, k)
    # if error code non-zero, a non-fatal error must have occurred
    # print warning, ignore this lambda value and return result
    # for all previous lambda values
    if (fit$jerr != 0) {
      errmsg <- jerr.multiview(fit$jerr, maxit, k)
      warning(errmsg$msg, call. = FALSE)
      k <- k - 1
      break
    }

    a0[k] <- fit$a0
    beta[, k] <- as.matrix(fit$beta)
    dev.ratio[k] <- fit$dev.ratio
    glambda[k] <- fit$lambda_scale
    
    # early stopping if dev.ratio almost 1 or no improvement
    if (k >= mnl && user_lambda == FALSE) {
      if (dev.ratio[k] > control$devmax) break
      else if (k > 1) {
        if (family$family == "gaussian") {
          if (dev.ratio[k] - dev.ratio[k-1] < control$fdev * dev.ratio[k])
            break
        } else if (family$family == "poisson") {
          if (dev.ratio[k] - dev.ratio[k - mnl + 1] <
                10 * control$fdev * dev.ratio[k])
            break
        } else if (dev.ratio[k] - dev.ratio[k-1] < control$fdev) break
      }
    }
  } # end of for(k in 1:nlam)
  
  if (trace.it == 1) {
    utils::setTxtProgressBar(pb, nlambda)
    cat("", fill = TRUE)
  }

  # truncate a0, beta, dev.ratio, lambda if necessary
  if (k < nlambda) {
    indices <- 1:k
    a0 <- a0[indices]
    beta <- beta[, indices, drop = FALSE]
    dev.ratio <- dev.ratio[indices]
    lambda <- lambda[indices]
    glambda <- glambda[indices]
  }

  ## So far glambda has merely been the scaling factor. Now fix it
  ## to reflect what it actually should be.
  glambda <- glambda * lambda
  
  # return coefficients to original scale (because of x standardization)
  beta <- beta / xs
  a0 <- a0 - colSums(beta * xm)

  # output
  stepnames <- paste0("s", 0:(length(lambda) - 1))
  out <- list()
  out$a0 <- a0
  names(out$a0) <- stepnames
  out$beta <- Matrix::Matrix(beta, sparse = TRUE,
                             dimnames = list(vnames, stepnames))
  out$df <- colSums(abs(beta) > 0)
  out$dim <- dim(beta)
  ## HERE IS A KEY SECTION OF CODE
  ## We always stick with the glmnet lambdas so as to concur with the case
  ## rho == 0. When rho == 0, we just call glmnet and the lambdas returned are
  ## the glmnet lambdas. They have the property of idempotence: call glmnet again with
  ## the returned lambda sequence returned and you get the same results.
  ## For rho > 0, this idempotence does not hold, which can be disconcerting!!
  ## This is because multiview.fit scales the lambda before calling glmnet:::elnet.
  ## So the lambda in the object should always be what glmnet routines were called with
  ## so as to use glmnet prediction methods etc. But we also store the unmodified lambdas
  ## as mvlambda, so that using mvlambda in the call will always guarantee idempotence
  ##
  out$lambda <- glambda
  out$mvlambda <- lambda
  ##
  ##
  out$dev.ratio <- dev.ratio
  out$nulldev <- start_val$nulldev
  out$npasses <- fit$npasses
  out$jerr <- fit$jerr
  out$offset <- is.offset
  ## out$call <- this.call
  out$family <- family
  out$nobs <- nobs
  ## ## We also need the standard deviations
  ## out$xsd <- xsd
  #class(out) <- c("multiview", "glmnetfit", "glmnet")
  class(out) <- c("glmnetfit", "glmnet")

  return(out)
}

#' Fit a GLM with elastic net regularization for a single value of
#' lambda
#'
#' Fit a generalized linear model via penalized maximum likelihood for
#' a single value of lambda. Can deal with any GLM family.
#'
#' WARNING: Users should not call `multiview.fit`
#' directly. Higher-level functions in this package call
#' `multiview.fit` as a subroutine. If a warm start object is
#' provided, some of the other arguments in the function may be
#' overriden.
#'
#' `multiview.fit` solves the elastic net problem for a _single,
#' user-specified_ value of lambda. `multiview.fit` works for any GLM
#' family. It solves the problem using iteratively reweighted least
#' squares (IRLS). For each IRLS iteration, `multiview.fit` makes a
#' quadratic (Newton) approximation of the log-likelihood, then calls
#' `elnet.fit` to minimize the resulting approximation.
#'
#' In terms of standardization: `multiview.fit` does not standardize
#' `x` and `weights`. `penalty.factor` is standardized so that to sum
#' to `nvars`.
#'
#' @inheritParams multiview.path
#' @param x the column-binded entries of `x_list`
#' @param lambda A single value for the `lambda` hyperparameter.
#' @param maxit Maximum number of passes over the data; default is
#'   `10^5`.  (If a warm start object is provided, the number of
#'   passes the warm start object performed is included.)
#' @param warm Either a `multiview` object or a list (with names
#'   `beta` and `a0` containing coefficients and intercept
#'   respectively) which can be used as a warm start. Default is
#'   `NULL`, indicating no warm start.  For internal use only.
#' @param from.multiview.path Was `multiview.fit()` called from
#'   `multiview.path()`?  Default is `FALSE`.This has implications for
#'   computation of the penalty factors.
#' @param save.fit Return the warm start object? Default is `FALSE`.
#' @param trace.it Controls how much information is printed to
#'   screen. If `trace.it = 2`, some information about the fitting
#'   procedure is printed to the console as the model is being
#'   fitted. Default is `trace.it = 0` (no information
#'   printed). (`trace.it = 1` not used for compatibility with
#'   `multiview.path`.)
#'
#' @return An object with class `"multiview"`. The list
#' returned contains more keys than that of a `"multiview"` object.
#' \item{a0}{Intercept value.}
#' \item{beta}{A `nvars` by `1` matrix of coefficients, stored in sparse matrix
#' format.}
#' \item{df}{The number of nonzero coefficients.}
#' \item{dim}{Dimension of coefficient matrix.}
#' \item{lambda}{Lambda value used.}
#' \item{lambda_scale}{The multiview lambda scale factor}
#' \item{dev.ratio}{The fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
#' for the saturated model (a model with a free parameter per observation).
#' Hence dev.ratio=1-dev/nulldev.}
#' \item{nulldev}{Null deviance (per observation). This is defined to be
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{A logical variable indicating whether an offset was included
#' in the model.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#' \item{warm_fit}{If `save.fit = TRUE`, output of C++ routine, used for
#' warm starts. For internal use only.}
#' \item{family}{Family used for the model.}
#' \item{converged}{A logical variable: was the algorithm judged to have
#' converged?}
#' \item{boundary}{A logical variable: is the fitted value on the boundary of
#' the attainable values?}
#' \item{obj_function}{Objective function value at the solution.}
#'
#' @importFrom stats gaussian
#' @importFrom utils combn
#' 
multiview.fit <- function(x_list, x, y, rho, weights, lambda, alpha = 1.0,
                          offset = rep(0, nobs), family = gaussian(),
                          intercept = TRUE, thresh = 1e-7, maxit = 100000,
                          penalty.factor = rep(1.0, nvars), exclude = c(), lower.limits = -Inf,
                          upper.limits = Inf, warm = NULL, from.multiview.path = FALSE,
                          save.fit = FALSE, trace.it = 0, user_lambda = FALSE) {
  this.call <- match.call()
  control <- multiview.control()

  nviews  <- length(x_list)
  p_x  <- lapply(x_list, ncol)
  ends  <- cumsum(p_x)
  starts  <- c(1, ends[-nviews] + 1)

  ### Prepare all the generic arguments
  nobs <- nrow(x)
  nvars <- ncol(x)
  is.offset <- !(missing(offset))
  if (is.offset == FALSE) {
    offset <- as.double(y * 0) #keeps the shape of y
  } else if (is.null(offset)) {
    offset <- rep(0, nobs)
    is.offset = FALSE
  }
 
  x_list <- lapply(split(seq_len(nvars), rep(seq_along(p_x), p_x)), function(ind) x[, ind])
  
  beta_indices <- mapply(seq.int, starts, ends, SIMPLIFY = FALSE)
  pairs <- apply(utils::combn(nviews, 2), 2, identity, simplify = FALSE)
  view_components <- lapply(pairs,
                            function(pair) {
                              i  <- pair[1L]; j  <- pair[2L];
                              list(index = list(beta_indices[[i]], beta_indices[[j]]),
                                   x = list(x_list[[i]], x_list[[j]]))
                            })

 # add xm and xs attributes if they are missing for sparse x
  # glmnet.fit assumes that x is already standardized. Any standardization
  # the user wants should be done beforehand.

  if (inherits(x, "sparseMatrix")) {
    if ("xm" %in% names(attributes(x)) == FALSE)
      attr(x, "xm") <- rep(0.0, times = nvars)
    if ("xs" %in% names(attributes(x)) == FALSE)
      attr(x, "xs") <- rep(1.0, times = nvars)
  }

  # if calling from glmnet.path(), we do not need to check on exclude
  # and penalty.factor arguments as they have been prepared by glmnet.path()
  if (!from.multiview.path) {
    # check and standardize penalty factors (to sum to nvars)
    if(any(penalty.factor == Inf)) {
      exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
      exclude = sort(unique(exclude))
    }
    if(length(exclude) > 0) {
      jd = match(exclude, seq(nvars), 0)
      if(!all(jd > 0)) stop ("Some excluded variables out of range")
      penalty.factor[jd] = 1 # ow can change lambda sequence
    }
    vp = pmax(0, penalty.factor)
    vp = as.double(vp * nvars / sum(vp))
  } else {
    vp <- as.double(penalty.factor)
  }

  ### check on limits
  ## lower.limits[lower.limits == -Inf] = -control$big
  ## upper.limits[upper.limits == Inf] = control$big
  ## if (length(lower.limits) < nvars)
  ##   lower.limits = rep(lower.limits, nvars) else
  ##                                             lower.limits = lower.limits[seq(nvars)]
  ## if (length(upper.limits) < nvars)
  ##   upper.limits = rep(upper.limits, nvars) else
  ##                                             upper.limits = upper.limits[seq(nvars)]
  ### end check on limits
  ### end preparation of generic arguments

  # get the relevant family functions
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv))
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE)
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if (is.null(x))
                                        if.null
  else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)

  # computation of null deviance (get mu in the process)
  if (is.null(warm)) {
    start_val <- get_start(x, y, weights, family, intercept, is.offset,
                           offset, exclude, vp, alpha)
    nulldev <- start_val$nulldev
    mu <- start_val$mu
    fit <- NULL
    coefold <- rep(0, nvars)   # initial coefs = 0
    eta <- family$linkfun(mu)
    intold <- (eta - offset)[1]
  } else {
    if ("glmnetfit" %in% class(warm)) {
      if (!inherits(warm$warm_fit, "warmfit")) stop("Invalid warm start object")
      fit <- warm
      nulldev <- fit$nulldev
      coefold <- fit$warm_fit$a   # prev value for coefficients
      intold <- fit$warm_fit$aint    # prev value for intercept
      eta <- get_eta(x, coefold, intold)
      mu <- linkinv(eta <- eta + offset)
    } else if ("list" %in% class(warm) && "a0" %in% names(warm) &&
                 "beta" %in% names(warm)) {
      nulldev <- get_start(x, y, weights, family, intercept, is.offset,
                           offset, exclude, vp, alpha)$nulldev
      fit <- warm
      coefold <- fit$beta   # prev value for coefficients
      intold <- fit$a0    # prev value for intercept
      eta <- get_eta(x, coefold, intold)
      mu <- linkinv(eta <- eta + offset)
    } else {
      stop("Invalid warm start object")
    }
  }

  if (!(validmu(mu) && valideta(eta)))
    stop("cannot find valid starting values: please specify some",
         call. = FALSE)

  start <- NULL     # current value for coefficients
  start_int <- NULL # current value for intercept
  obj_val_old <- obj_function(y, mu, weights, family, lambda, alpha, coefold, vp, view_components, rho)
  if (trace.it == 2) {
    cat("Warm Start Objective:", obj_val_old, fill = TRUE)
  }

  ## precompute the fixed offset vector for the larger problem
  g_offset <- c(offset, rep(offset, length(pairs)))  #NOTE!!


  conv <- FALSE      # converged?

  sum_weights <- sum(weights)
  # IRLS loop
  for (iter in 1L:control$mxitnr) {
    # some checks for NAs/zeros
    xx  <- x
    varmu <- variance(mu)
    if (anyNA(varmu)) stop("NAs in V(mu)")
    if (any(varmu == 0)) stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")

    # compute working response and weights
    zz <- (eta - offset) + (y - mu)/mu.eta.val
    w <- (weights * mu.eta.val^2)/variance(mu)
    # have to update the weighted residual in our fit object
    # (in theory g and iy should be updated too, but we actually recompute g
    # and iy anyway in wls.f)
    ## if (!is.null(fit)) {
    ##   fit$warm_fit$r <- w * (zz - eta + offset)
    ## }

    w_sum <- sum(w)
    w_std <- w / sum(w)

    mzz <- sum(w * zz) / w_sum
    zzc <- zz - mzz

    mx <- apply(xx, 2, function(x) sum(w_std * x))

    if (!inherits(xx, "sparseMatrix")) {
      xx  <- sweep(xx, 2L, mx, check.margin = FALSE)
    } else {
      attr(xx, "xm") <- mx
      attr(xx, "xs") <- rep(1.0, times = nvars)
    }

    nx_list <- lapply(x_list, function(mat) {
      column_means <- apply(mat, 2, function(column) sum(w_std * column))
      sweep(mat, 2L, column_means, check.margin = FALSE)
    })
    features <- xx
    target <- zzc
    rows <- lapply(pairs, make_row, x_list = nx_list, p_x = p_x, rho = rho )
    features <- do.call(rbind, c(list(features), rows))
    target <- c(target, rep(0, length(pairs) * nobs))
    w <- c(w, rep(weights, length(pairs)))  #NOTE!!
    w_sum <- sum(w)
    w_std <- w / w_sum
    if (!is.null(fit)) {
      ## features and w_std below will have the right dimension from previous iteration!
      g_offset  <- c(offset, rep(offset, length(pairs)))
      g_eta <- get_eta(features, fit$warm$a, 0) ## intercept is zero for larger fit!
      fit$warm_fit$r <- w_std * (target - g_eta + g_offset)
    }

    #cat(sprintf("SW is %f\n", w_sum))
    ## NOTE: sum(weights) below takes care of glmnet parameterization lambda -> lambda * n!
    if (user_lambda) {
      lambda2 <- lambda
    } else {
      lambda2  <- lambda * sum_weights / w_sum
    }
    #print(w)
    #cat(sprintf("Sum wt: %f Our Lambda%f\n", w_sum, lambda2))

    ## fit <- glmnet:::glmnet.fit(features, target, lambda = lambda2, weights = w_std,
    ##                            intercept = FALSE, from.glmnet.path = TRUE, save.fit = TRUE, trace.it = 2)

    ## TODO: We don't use warm starts yet!
    ## To do that, we need to recompute the residuals else the fit fails
    ## Which is why `warm = fit` is commented out below.
    fit <- elnet.fit(x = features, y = target, weights = w_std,
                     lambda = lambda2, alpha = alpha,
                     exclude = exclude,
                     intercept = FALSE, from.glmnet.fit = TRUE, save.fit = TRUE,
                     thresh = thresh, maxit = maxit,
                     upper.limits = upper.limits, penalty.factor = vp, warm = fit
                     )

    if (fit$jerr != 0) return(list(jerr = fit$jerr))

    # update coefficients, eta, mu and obj_val
    start <- fit$warm_fit$a
    #start_int <- fit$warm_fit$aint
    start_int <- mzz - sum(mx * start)
    eta <- get_eta(x, start, start_int)
    mu <- linkinv(eta <- eta + offset)
    obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp, view_components, rho)
    if (trace.it == 2) cat("Iteration", iter, "Objective:", obj_val, fill = TRUE)

    boundary <- FALSE
    halved <- FALSE  # did we have to halve the step size?
    # if objective function is not finite, keep halving the stepsize until it is finite
    # for the halving step, we probably have to adjust fit$g as well?
    if (!is.finite(obj_val) || obj_val > control$big) {
      warning("Infinite objective function!", call. = FALSE)
      if (is.null(coefold) || is.null(intold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      warning("step size truncated due to divergence", call. = FALSE)
      ii <- 1
      while (!is.finite(obj_val) || obj_val > control$big) {
        if (ii > control$mxitnr)
          stop("inner loop 1; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        start_int <- (start_int + intold)/2
        eta <- get_eta(x, start, start_int)
        mu <- linkinv(eta <- eta + offset)
        obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp, view_components, rho)
        if (trace.it == 2) cat("Iteration", iter, " Halved step 1, Objective:",
                               obj_val, fill = TRUE)
      }
      boundary <- TRUE
      halved <- TRUE
    }
    # if some of the new eta or mu are invalid, keep halving stepsize until valid
    if (!(valideta(eta) && validmu(mu))) {
      warning("Invalid eta/mu!", call. = FALSE)
      if (is.null(coefold) || is.null(intold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      warning("step size truncated: out of bounds", call. = FALSE)
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > control$mxitnr)
          stop("inner loop 2; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        start_int <- (start_int + intold)/2
        eta <- get_eta(x, start, start_int)
        mu <- linkinv(eta <- eta + offset)
      }
      boundary <- TRUE
      halved <- TRUE
      obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp, view_components, rho)
      if (trace.it == 2) cat("Iteration", iter, " Halved step 2, Objective:", obj_val, fill = TRUE)
    }
    # extra halving step if objective function value actually increased
    if (obj_val > obj_val_old + 1e-7) {
      ii <- 1
      while (obj_val > obj_val_old + 1e-7) {
        ##cat(sprintf("Iter: %d, Diff: %10.f\n", ii, obj_val - obj_val_old))
        if (ii > control$mxitnr)
          stop("inner loop 3; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        start_int <- (start_int + intold)/2
        eta <- get_eta(x, start, start_int)
        mu <- linkinv(eta <- eta + offset)
        obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp, view_components, rho)
        if (trace.it == 2) cat(sprintf("Iteration %d, Halved step 3, Objective: %.10f\n", iter, obj_val))
      }
      halved <- TRUE
    }

    # if we did any halving, we have to update the coefficients, intercept
    # and weighted residual in the warm_fit object
    if (halved) {
      fit$warm_fit$a <- start
      fit$warm_fit$aint <- start_int
      g_eta <- get_eta(features, start, 0) ## intercept is zero for larger fit!
      fit$warm_fit$r <- w_std * (target - g_eta) + g_offset
    }

    # test for convergence
    if (abs(obj_val - obj_val_old)/(0.1 + abs(obj_val)) < control$epsnr) {
      conv <- TRUE
      break
    }
    else {
      coefold <- start
      intold <- start_int
      obj_val_old <- obj_val
    }
  }
  # end of IRLS loop

  ## Fix up a0 for coeff!
  fit$a0  <- start_int
  ## The scale below is used to determine the actual lambda seq for multiview
  if (user_lambda) {
    fit$lambda_scale <- 1
  } else {
    fit$lambda_scale <- sum_weights / w_sum
  }
  # checks on convergence and fitted values
  if (!conv)
    warning("glmnet.fit: algorithm did not converge", call. = FALSE)
  if (boundary)
    warning("glmnet.fit: algorithm stopped at boundary value", call. = FALSE)

  # some extra warnings, printed only if trace.it == 2
  if (trace.it == 2) {
    eps <- 10 * .Machine$double.eps
    if ((family$family == "binomial") && (any(mu > 1 - eps) || any(mu < eps)))
      warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
              call. = FALSE)
    if ((family$family == "poisson") && (any(mu < eps)))
      warning("glm.fit: fitted rates numerically 0 occurred",
              call. = FALSE)
  }

  # prepare output object
  if (save.fit == FALSE) {
    fit$warm_fit <- NULL
  }

  # overwrite values from elnet.fit object
  fit$call <- this.call
  fit$offset <- is.offset
  fit$nulldev <- nulldev
  fit$dev.ratio <- 1 - dev_function(y, mu, weights, family) / nulldev
  ##fit$dev.ratio <- 1 - dev_function(y, mu, w_std, family) / nulldev

  # add new key-value pairs to list
  fit$family <- family
  fit$converged <- conv
  fit$boundary <- boundary
  fit$obj_function <- obj_val

  class(fit) <- c("multiview", "glmnetfit", "glmnet")
  fit
}


jerr.multiview <- function (n, maxit, k = NULL) {
    if (n == 0) {
        list(n = 0, fatal = FALSE, msg = "")
    } else if (n > 0) {
        # fatal error
        fatal <- TRUE
        msg <- ifelse(n < 7777,
                      "Memory allocation error; contact package maintainer",
                      "Unknown error")
    } else {
        # non-fatal error
        fatal <- FALSE
        msg <- paste("Convergence for ", k, "th lambda value not reached after maxit=",
                     maxit, " iterations; solutions for larger lambdas returned",
                     sep = "")
    }
    list(n = n, fatal = fatal, msg = msg)
}

#' Solve weighted least squares (WLS) problem for a single lambda value
#'
#' Solves the weighted least squares (WLS) problem for a single lambda value. Internal
#' function that users should not call directly.
#'
#' WARNING: Users should not call \code{elnet.fit} directly. Higher-level functions
#' in this package call \code{elnet.fit} as a subroutine. If a warm start object
#' is provided, some of the other arguments in the function may be overriden.
#'
#' \code{elnet.fit} is essentially a wrapper around a C++ subroutine which
#' minimizes
#'
#' \deqn{1/2 \sum w_i (y_i - X_i^T \beta)^2 + \sum \lambda \gamma_j
#' [(1-\alpha)/2 \beta^2+\alpha|\beta|],}
#'
#' over \eqn{\beta}, where \eqn{\gamma_j} is the relative penalty factor on the
#' jth variable. If \code{intercept = TRUE}, then the term in the first sum is
#' \eqn{w_i (y_i - \beta_0 - X_i^T \beta)^2}, and we are minimizing over both
#' \eqn{\beta_0} and \eqn{\beta}.
#'
#' None of the inputs are standardized except for \code{penalty.factor}, which
#' is standardized so that they sum up to \code{nvars}.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed that any standardization needed
#' has already been done.
#' @param y Quantitative response variable.
#' @param weights Observation weights. \code{elnet.fit} does NOT standardize
#' these weights.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' The penalty is defined as \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
#' \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)?
#' @param thresh Convergence threshold for coordinate descent. Each inner
#' coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than thresh times the null deviance.
#' Default value is \code{1e-7}.
#' @param maxit Maximum number of passes over the data; default is \code{10^5}.
#' (If a warm start object is provided, the number of passes the warm start object
#' performed is included.)
#' @param penalty.factor Separate penalty factors can be applied to each
#' coefficient. This is a number that multiplies \code{lambda} to allow differential
#' shrinkage. Can be 0 for some variables, which implies no shrinkage, and that
#' variable is always included in the model. Default is 1 for all variables (and
#' implicitly infinity for variables listed in exclude). Note: the penalty
#' factors are internally rescaled to sum to \code{nvars}.
#' @param exclude Indices of variables to be excluded from the model. Default is
#' none. Equivalent to an infinite penalty factor.
#' @param lower.limits Vector of lower limits for each coefficient; default
#' \code{-Inf}. Each of these must be non-positive. Can be presented as a single
#' value (which will then be replicated), else a vector of length \code{nvars}.
#' @param upper.limits Vector of upper limits for each coefficient; default
#' \code{Inf}. See \code{lower.limits}.
#' @param warm Either a \code{glmnetfit} object or a list (with names \code{beta}
#' and \code{a0} containing coefficients and intercept respectively) which can
#' be used as a warm start. Default is \code{NULL}, indicating no warm start.
#' For internal use only.
#' @param from.glmnet.fit Was \code{elnet.fit()} called from \code{glmnet.fit()}?
#' Default is FALSE.This has implications for computation of the penalty factors.
#' @param save.fit Return the warm start object? Default is FALSE.
#'
#' @return An object with class "glmnetfit" and "glmnet". The list returned has
#' the same keys as that of a \code{glmnet} object, except that it might have an
#' additional \code{warm_fit} key.
#' \item{a0}{Intercept value.}
#' \item{beta}{A \code{nvars x 1} matrix of coefficients, stored in sparse matrix
#' format.}
#' \item{df}{The number of nonzero coefficients.}
#' \item{dim}{Dimension of coefficient matrix.}
#' \item{lambda}{Lambda value used.}
#' \item{dev.ratio}{The fraction of (null) deviance explained. The deviance
#' calculations incorporate weights if present in the model. The deviance is
#' defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
#' for the saturated model (a model with a free parameter per observation).
#' Hence dev.ratio=1-dev/nulldev.}
#' \item{nulldev}{Null deviance (per observation). This is defined to be
#' 2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
#' \item{npasses}{Total passes over the data.}
#' \item{jerr}{Error flag, for warnings and errors (largely for internal
#' debugging).}
#' \item{offset}{Always FALSE, since offsets do not appear in the WLS problem.
#' Included for compability with glmnet output.}
#' \item{call}{The call that produced this object.}
#' \item{nobs}{Number of observations.}
#' \item{warm_fit}{If \code{save.fit=TRUE}, output of C++ routine, used for
#' warm starts. For internal use only.}
#'
elnet.fit <- function(x, y, weights, lambda, alpha = 1.0, intercept = TRUE,
                      thresh = 1e-7, maxit = 100000,
                      penalty.factor = rep(1.0, nvars), exclude = c(),
                      lower.limits = -Inf, upper.limits = Inf, warm = NULL,
                      from.glmnet.fit = FALSE, save.fit = FALSE) {
    this.call <- match.call()
    internal.parms <- multiview.control()

    # compute null deviance
    ybar <- weighted.mean(y, weights)
    nulldev <- sum(weights * (y - ybar)^2)

    # if class "glmnetfit" warmstart object provided, pull whatever we want out of it
    # else, prepare arguments, then check if coefs provided as warmstart
    # (if only coefs are given as warmstart, we prepare the other arguments
    # as if no warmstart was provided)
    if (!is.null(warm) && "glmnetfit" %in% class(warm)) {
        warm <- warm$warm_fit
        if (!inherits(warm, "warmfit")) stop("Invalid warm start object")

        a <- warm$a
        aint <- warm$aint
        alm0 <- warm$almc
        cl <- warm$cl
        g <- warm$g
        ia <- warm$ia
        iy <- warm$iy
        iz <- warm$iz
        ju <- warm$ju
        m <- warm$m
        mm <- warm$mm
        nino <- warm$nino
        nobs <- warm$no
        nvars <- warm$ni
        nlp <- warm$nlp
        nx <- warm$nx
        r <- warm$r
        rsqc <- warm$rsqc
        xv <- warm$xv
        vp <- warm$vp
    } else {
        nobs <- as.integer(nrow(x))
        nvars <- as.integer(ncol(x))

        # if calling from glmnet.fit(), we do not need to check on exclude
        # and penalty.factor arguments as they have been prepared by glmnet.fit()
        # Also exclude will include variance 0 columns
        if (!from.glmnet.fit) {
            # check and standardize penalty factors (to sum to nvars)
            if(any(penalty.factor == Inf)) {
                exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
                exclude = sort(unique(exclude))
            }
            if(length(exclude) > 0) {
                jd = match(exclude, seq(nvars), 0)
                if(!all(jd > 0)) stop ("Some excluded variables out of range")
                penalty.factor[jd] = 1 # ow can change lambda sequence
            }
            vp = pmax(0, penalty.factor)
            vp = as.double(vp * nvars / sum(vp))
        } else {
            vp <- as.double(penalty.factor)
        }
        # compute ju
        # assume that there are no constant variables
        ju <- rep(1, nvars)
        ju[exclude] <- 0
        ju <- as.integer(ju)

        # compute cl from lower.limits and upper.limits
        lower.limits[lower.limits == -Inf] <- -internal.parms$big
        upper.limits[upper.limits == Inf] <- internal.parms$big
        if (length(lower.limits) < nvars)
            lower.limits = rep(lower.limits, nvars) else
                lower.limits = lower.limits[seq(nvars)]
        if (length(upper.limits) < nvars)
            upper.limits = rep(upper.limits, nvars) else
                upper.limits = upper.limits[seq(nvars)]
        cl <- rbind(lower.limits, upper.limits)
        storage.mode(cl) = "double"

        nx <- as.integer(nvars)

        a <- double(nvars)
        aint <- double(1)
        alm0 <- double(1)
        g <- double(nvars)
        ia <- integer(nx)
        iy <- integer(nvars)
        iz <- integer(1)
        m <- as.integer(1)
        mm <- integer(nvars)
        nino <- integer(1)
        nlp <- integer(1)
        r <- weights * y
        rsqc <- double(1)
        xv <- double(nvars)

        # check if coefs were provided as warmstart: if so, use them
        if (!is.null(warm)) {
            if ("list" %in% class(warm) && "a0" %in% names(warm) &&
                "beta" %in% names(warm)) {
                a <- as.double(warm$beta)
                aint <- as.double(warm$a0)
                mu <- drop(x %*% a + aint)
                r <- weights * (y - mu)
                rsqc <- 1 - sum(weights * (y - mu)^2) / nulldev
            } else {
                stop("Invalid warm start object")
            }
        }
    }

    # for the parameters here, we are overriding the values provided by the
    # warmstart object
    alpha <- as.double(alpha)
    almc <- as.double(lambda)
    intr <- as.integer(intercept)
    jerr <- integer(1)
    maxit <- as.integer(maxit)
    thr <- as.double(thresh)
    v <- as.double(weights)

    a.new <- a
    a.new[0] <- a.new[0] # induce a copy
  
    # take out components of x and run C++ subroutine
    if (inherits(x, "sparseMatrix")) {
        xm <- as.double(attr(x, "xm"))
        xs <- as.double(attr(x, "xs"))
        wls_fit <- spwls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
                             x=x,xm=xm,xs=xs,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
                             maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
                             nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
    } else {
        wls_fit <- wls_exp(alm0=alm0,almc=almc,alpha=alpha,m=m,no=nobs,ni=nvars,
                           x=x,r=r,xv=xv,v=v,intr=intr,ju=ju,vp=vp,cl=cl,nx=nx,thr=thr,
                           maxit=maxit,a=a.new,aint=aint,g=g,ia=ia,iy=iy,iz=iz,mm=mm,
                           nino=nino,rsqc=rsqc,nlp=nlp,jerr=jerr)
    }

    # if error code > 0, fatal error occurred: stop immediately
    # if error code < 0, non-fatal error occurred: return error code
    if (wls_fit$jerr > 0) {
        errmsg <- jerr.glmnetfit(wls_fit$jerr, maxit)
        stop(errmsg$msg, call. = FALSE)
    } else if (wls_fit$jerr < 0)
        return(list(jerr = wls_fit$jerr))
    warm_fit <- wls_fit[c("almc", "r", "xv", "ju", "vp", "cl", "nx",
                          "a", "aint", "g", "ia", "iy", "iz", "mm", "nino",
                          "rsqc", "nlp")]
    warm_fit[['m']] <- m
    warm_fit[['no']] <- nobs
    warm_fit[['ni']] <- nvars
    class(warm_fit) <- "warmfit"

    beta <- Matrix::Matrix(wls_fit$a, sparse = TRUE)

    out <- list(a0 = wls_fit$aint, beta = beta, df = sum(abs(beta) > 0),
                dim = dim(beta), lambda = lambda, dev.ratio = wls_fit$rsqc,
                nulldev = nulldev, npasses = wls_fit$nlp, jerr = wls_fit$jerr,
                offset = FALSE, call = this.call, nobs = nobs, warm_fit = warm_fit)
    if (save.fit == FALSE) {
        out$warm_fit <- NULL
    }
    class(out) <- c("glmnetfit", "glmnet")
    out
}

