#' Perform cooperative learning using the direct algorithm for
#' two or more views.
#'
#' `multiview` uses [glmnet::glmnet()] to do most of its work and
#' therefore takes many of the same parameters, but an intercept is
#' always included and several other parameters do not
#' apply. Such inapplicable arguments are overridden and warnings issued.
#' 
#' The current code can be slow for "large" data sets, e.g. when the
#' number of features is larger than 1000.  It can be helpful to see
#' the progress of multiview as it runs; to do this, set trace.it = 1
#' in the call to multiview or cv.multiview.  With this, multiview
#' prints out its progress along the way.  One can also pre-filter the
#' features to a smaller set, using the exclude option, with a filter
#' function.
#' 
#' If there are missing values in the feature matrices: we recommend
#' that you center the columns of each feature matrix, and then fill
#' in the missing values with 0.
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
#' @param x_list a list of `x` matrices with same number of rows
#'   `nobs`
#' @param y the quantitative response with length equal to `nobs`, the
#'   (same) number of rows in each `x` matrix
#' @param rho the weight on the agreement penalty, default 0. `rho=0`
#'   is a form of early fusion, and `rho=1` is a form of late fusion.
#'   We recommend trying a few values of `rho` including 0, 0.1, 0.25,
#'   0.5, and 1 first; sometimes `rho` larger than 1 can also be
#'   helpful.
#' @param family A description of the error distribution and link
#'   function to be used in the model. This is the result of a call to
#'   a family function. Default is [stats::gaussian]. (See
#'   [stats::family] for details on family functions.)
#' @param weights observation weights. Can be total counts if
#'   responses are proportion matrices. Default is 1 for each
#'   observation
#' @param offset A vector of length \code{nobs} that is included in
#'   the linear predictor (a \code{nobs x nc} matrix for the
#'   \code{"multinomial"} family).  Useful for the \code{"poisson"}
#'   family (e.g. log of exposure time), or for refining a model by
#'   starting at a current fit. Default is \code{NULL}. If supplied,
#'   then values must also be supplied to the \code{predict} function.
#' @param alpha The elasticnet mixing parameter, with
#'   \eqn{0\le\alpha\le 1}.  The penalty is defined as
#'   \deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
#'   \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge
#'   penalty.
#' @param nlambda The number of \code{lambda} values - default is 100.
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a
#'   fraction of \code{lambda.max}, the (data derived) entry value
#'   (i.e. the smallest value for which all coefficients are
#'   zero). The default depends on the sample size \code{nobs}
#'   relative to the number of variables \code{nvars}. If \code{nobs >
#'   nvars}, the default is \code{0.0001}, close to zero.  If
#'   \code{nobs < nvars}, the default is \code{0.01}.  A very small
#'   value of \code{lambda.min.ratio} will lead to a saturated fit in
#'   the \code{nobs < nvars} case. This is undefined for
#'   \code{"binomial"} and \code{"multinomial"} models, and
#'   \code{glmnet} will exit gracefully when the percentage deviance
#'   explained is almost 1.
#' @param lambda A user supplied \code{lambda} sequence. Typical usage
#'   is to have the program compute its own \code{lambda} sequence
#'   based on \code{nlambda} and \code{lambda.min.ratio}. Supplying a
#'   value of \code{lambda} overrides this. WARNING: use with
#'   care. Avoid supplying a single value for \code{lambda} (for
#'   predictions after CV use \code{predict()} instead).  Supply
#'   instead a decreasing sequence of \code{lambda}
#'   values. \code{glmnet} relies on its warms starts for speed, and
#'   its often faster to fit a whole path than compute a single fit.
#' @param standardize Logical flag for x variable standardization,
#'   prior to fitting the model sequence. The coefficients are always
#'   returned on the original scale. Default is
#'   \code{standardize=TRUE}.  If variables are in the same units
#'   already, you might not wish to standardize. See details below for
#'   y standardization with \code{family="gaussian"}.
#' @param intercept Should intercept(s) be fitted (default `TRUE`)
#' @param thresh Convergence threshold for coordinate descent. Each
#'   inner coordinate-descent loop continues until the maximum change
#'   in the objective after any coefficient update is less than
#'   \code{thresh} times the null deviance. Defaults value is
#'   \code{1E-7}.
#' @param penalty.factor Separate penalty factors can be applied to
#'   each coefficient. This is a number that multiplies \code{lambda}
#'   to allow differential shrinkage. Can be 0 for some variables,
#'   which implies no shrinkage, and that variable is always included
#'   in the model. Default is 1 for all variables (and implicitly
#'   infinity for variables listed in \code{exclude}). Note: the
#'   penalty factors are internally rescaled to sum to nvars, and the
#'   lambda sequence will reflect this change.
#' @param lower.limits Vector of lower limits for each coefficient;
#'   default \code{-Inf}. Each of these must be non-positive. Can be
#'   presented as a single value (which will then be replicated), else
#'   a vector of length \code{nvars}
#' @param upper.limits Vector of upper limits for each coefficient;
#'   default \code{Inf}. See \code{lower.limits}
#' @param maxit Maximum number of passes over the data for all lambda
#'   values; default is 10^5.
#' @param exclude Indices of variables to be excluded from the
#'   model. Default is none. Equivalent to an infinite penalty factor
#'   for the variables excluded (next item).  Users can supply instead
#'   an `exclude` function that generates the list of indices.  This
#'   function is most generally defined as `function(x_list, y, ...)`,
#'   and is called inside `multiview` to generate the indices for
#'   excluded variables.  The `...` argument is required, the others
#'   are optional.  This is useful for filtering wide data, and works
#'   correctly with `cv.multiview`. See the vignette 'Introduction'
#'   for examples.
#' @param lambda A user supplied `lambda` sequence, default
#'   `NULL`. Typical usage is to have the program compute its own
#'   `lambda` sequence. This sequence, in general, is different from
#'   that used in the [glmnet::glmnet()] call (named `lambda`)
#'   Supplying a value of `lambda` overrides this. WARNING: use with
#'   care. Avoid supplying a single value for `lambda` (for
#'   predictions after CV use [stats::predict()] instead.  Supply
#'   instead a decreasing sequence of `lambda` values as `multiview`
#'   relies on its warms starts for speed, and its often faster to fit
#'   a whole path than compute a single fit.
#' @param trace.it If \code{trace.it=1}, then a progress bar is
#'   displayed; useful for big models that take a long time to fit.
#' @return An object with S3 class `"multiview","*" `, where `"*"` is
#'   `"elnet"`, `"lognet"`, `"multnet"`, `"fishnet"` (poisson),
#'   `"coxnet"` or `"mrelnet"` for the various types of models.
#'   \item{call}{the call that produced this object}
#'   \item{a0}{Intercept sequence of length `length(lambda)`}
#'   \item{beta}{For `"elnet"`, `"lognet"`, `"fishnet"` and `"coxnet"`
#'   models, a `nvars x length(lambda)` matrix of coefficients, stored
#'   in sparse column format (`"CsparseMatrix"`). For `"multnet"` and
#'   `"mgaussian"`, a list of `nc` such matrices, one for each class.}
#'   \item{lambda}{The actual sequence of [glmnet::glmnet()] `lambda`
#'   values used. When `alpha=0`, the largest lambda reported does not
#'   quite give the zero coefficients reported (`lambda=inf` would in
#'   principle).  Instead, the largest `lambda` for `alpha=0.001` is
#'   used, and the sequence of `lambda` values is derived from this.}
#'   \item{lambda}{The sequence of lambda values} \item{mvlambda}{The
#'   corresponding sequence of multiview lambda values}
#'   \item{dev.ratio}{The fraction of (null) deviance explained (for
#'   `"elnet"`, this is the R-square). The deviance calculations
#'   incorporate weights if present in the model. The deviance is
#'   defined to be 2*(loglike_sat - loglike), where loglike_sat is the
#'   log-likelihood for the saturated model (a model with a free
#'   parameter per observation). Hence dev.ratio=1-dev/nulldev.}
#'   \item{nulldev}{Null deviance (per observation). This is defined
#'   to be 2*(loglike_sat -loglike(Null)); The NULL model refers to
#'   the intercept model, except for the Cox, where it is the 0
#'   model.} \item{df}{The number of nonzero coefficients for each
#'   value of `lambda`. For `"multnet"`, this is the number of
#'   variables with a nonzero coefficient for \emph{any} class.}
#'   \item{dfmat}{For `"multnet"` and `"mrelnet"` only. A matrix
#'   consisting of the number of nonzero coefficients per class}
#'   \item{dim}{dimension of coefficient matrix (ices)}
#'   \item{nobs}{number of observations} \item{npasses}{total passes
#'   over the data summed over all lambda values} \item{offset}{a
#'   logical variable indicating whether an offset was included in the
#'   model} \item{jerr}{error flag, for warnings and errors (largely
#'   for internal debugging).}
#'   
#' @seealso \code{print}, \code{coef}, \code{coef_ordered}, \code{predict}, and \code{plot}
#' methods for \code{"multiview"}, and the \code{"cv.multiview"} function.
#' 
#' @examples
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = matrix(rnorm(100 * 10), 100, 10)
#' y = rnorm(100)
#' fit1 = multiview(list(x=x,z=z), y, rho = 0)
#' print(fit1)
#' 
#' # extract coefficients at a single value of lambda
#' coef(fit1, s = 0.01) 
#' 
#' # extract ordered (standardized) coefficients at a single value of lambda
#' coef_ordered(fit1, s = 0.01) 
#' 
#' # make predictions
#' predict(fit1, newx = list(x[1:10, ],z[1:10, ]), s = c(0.01, 0.005))
#' 
#' # make a path plot of features for the fit
#' plot(fit1, label=TRUE)
#' 
#' # Binomial
#' by = sample(c(0,1), 100, replace = TRUE)
#' fit2 = multiview(list(x=x,z=z), by, family = binomial(), rho=0.5)
#' predict(fit2, newx = list(x[1:10, ],z[1:10, ]), s = c(0.01, 0.005), type="response")
#' coef_ordered(fit2, s = 0.01) 
#' plot(fit2, label=TRUE)
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = multiview(list(x=x,z=z), py, family = poisson(), rho=0.5)
#' predict(fit3, newx = list(x[1:10, ],z[1:10, ]), s = c(0.01, 0.005), type="response")
#' coef_ordered(fit3, s = 0.01) 
#' plot(fit3, label=TRUE)
#'
#' @importFrom glmnet glmnet
#' @importFrom stats family gaussian sd weighted.mean
#' @export
multiview <- function(x_list, y, rho = 0, family = gaussian(), weights = NULL, offset = NULL,
                      alpha = 1.0, nlambda = 100, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001),
                      lambda = NULL, standardize = TRUE, intercept = TRUE, thresh = 1e-7,
                      maxit = 100000, penalty.factor = rep(1.0, nvars), exclude = list(), lower.limits = -Inf,
                      upper.limits = Inf, trace.it = 0) {

  ## TODO: alpha = 1.0, stdize = TRUE, intercept = TRUE, exclude if list of ints have to be checked
  this.call <- match.call()

  ## TODO: CLEAN up
  mvlambda <- lambda
  is.offset <- !(is.null(offset))

  cox_case <- gaussian_case <- FALSE
  if (is.function(family)) {
    family <- family()
  }
  if (inherits(family, "family")) {
    gaussian_case <- (family$family == "gaussian" && family$link == "identity")
  } else if (is.character(family) && family == "cox") {
    cox_case <- TRUE
  } else {
    stop("multiview: family must be a family function or the string 'cox'")
  }

  ### Need to do this first so defaults in call can be satisfied

  if (rho < 0 || length(row) != 1L) {
    stop("multiview: rho must be a scalar >= 0!")
  }

  nviews  <- length(x_list)
  if (nviews < 1L) {
     stop("multiview: need at least one view!")
  }

  ## We cannot do this check because sometimes view_contribution calls multiview with some columns removed
  ## if (nviews == 1L && rho > 0) {
  ##   stop("multiview: a nonzero rho requires more than one view.")
  ## }

  dims  <- lapply(x_list, dim)
  if (any(sapply(dims, is.null))) {
    stop("multiview: x_list should be a list of matrices!")
  }

  nobs <- dims[[1L]][1L]
  if (any(sapply(x_list, nrow) != nobs)) {
    stop("multiview: all views must have the same number of rows!")
  }

  ## Give the x_list components names
  x_names  <- names(x_list)
  if (is.null(x_names) || (any(x_names == ""))) {
    new_names  <- sprintf("View%d", seq_len(nviews))
    x_names <- names(x_list)  <- new_names
  }

  p_x <- lapply(x_list, ncol)
  ## We also need colnames of the x matrices for decipherable output later, so we store them
  colnames_list <- lapply(seq_along(x_list), function(i) {
    cn <- colnames(x_list[[i]])
    if (is.null(cn)) sprintf("V%d", seq_len(p_x[[i]])) else cn
  })

  nvars <- Reduce(f = sum, p_x)

  ## ## Prepare to reuse glmnetFlex code
  ## x <- do.call(cbind, x_list)
  ## ## We need the std devs for other purposes, so we compute it
  ## ## TODO: Make this take weights in to account
  ## xsd <- apply(x, 2, sd)

  ## Handle exclude like in glmnet
  if (is.null(exclude) || length(exclude) == 0) {
    exclude <- integer(0)
  } else {
    if (is.function(exclude)) {
      exclude <- exclude(x_list = x_list, y = y, weights = weights)
    }
    exclude <- process_exclude(exclude, p_x)
  }
  
  if (rho > 0 && nviews > 1L) {
    ## Compute some things as we reuse them over and over again
    alpha <- as.double(alpha)
    
    # get feature variable names
    vnames <- unlist(mapply(paste, names(x_list), colnames_list, sep=":"))
    
    #if(is.null(vnames))vnames=paste("V",seq(nvars),sep="")
    
    # check weights
    if(is.null(weights)) {
      weights <- rep(1.0, nobs)
    } else {
      if (length(weights) != nobs)
        stop(paste("Number of elements in weights (",length(weights),
                   ") not equal to the number of rows of x (",nobs,")",sep=""))
      weights <- as.double(weights)
    }
    
    ## initialize from family function. Makes y a vector in case of binomial, and possibly changes weights
    ## Expects nobs to be defined, and creates n and mustart (neither used here)
    ## Some cases expect to see things, so we set it up just to make it work
    if (!is.character(family)) {
      etastart <- 0; mustart <- NULL; start <- NULL
      eval(family$initialize)
      
      ## Also make sure intercept is true
      if (!intercept) {
        ## BN TODO: Why not remove intercept as an option since we totally control it??
        warning("multiview: overriding intercept to TRUE for glm families")
        intercept <- TRUE
      }
    } else {
      ## We ignore intercept for Cox
      intercept <- FALSE
    }
    
    if (!is.offset) {
      #offset <- as.double(y * 0) #keeps the shape of y
      #offset <- rep(0, times = nrow(y))  # for cox
      offset <- numeric(nobs)
    }
    
    # infinite penalty factor vars are excluded
    exclude <- c(exclude, seq(nvars)[penalty.factor == Inf])
    exclude <- sort(unique(exclude))
    
    ## Compute weighted mean and variance of columns of x, sensitive to sparse matrix
    ## needed to detect constant columns below, and later if standarization
    ##meansd <- weighted_mean(x, weights)

    mv_meansd <- lapply(x_list, weighted_mean_sd, weights = weights)
    col_means_list <- lapply(mv_meansd, function(x) x$mean)
    col_means <- unlist(col_means_list)
    col_sd_list <- lapply(mv_meansd, function(x) x$sd)
    col_sd <- unlist(col_sd_list)
    
    ## look for constant variables, and if any, then add to exclude
    const_vars <- col_sd == 0
    nzvar <- setdiff(which(!const_vars), exclude)
    # if all the non-excluded variables have zero variance, throw error
    if (length(nzvar) == 0) stop("All used predictors have zero variance")
    
    ## if any constant vars, add to exclude
    if(any(const_vars)) {
      exclude <- sort(unique(c(which(const_vars),exclude)))
      col_sd[const_vars] <- 1.0 ## we divide later, and do not want bad numbers
      col_sd_list <- lapply(col_sd_list, function(x) ifelse(x == 0, 1.0, x)) # shadowing col_sd!
    }
    if(length(exclude) > 0) {
      jd <- match(exclude, seq(nvars), 0)
      if(!all(jd > 0)) stop ("Some excluded variables out of range")
      penalty.factor[jd] <- 1 # ow can change lambda sequence
    }
    # check and standardize penalty factors (to sum to nvars)
    vp <- pmax(0, penalty.factor)
    if (max(vp) <= 0) stop("All penalty factors are <= 0")
    vp <- as.double(vp * nvars / sum(vp))
    
    ### check on limits
    control <- multiview.control()
    if (thresh >= control$epsnr)
      warning("thresh should be smaller than multiview.control()$epsnr",
              call. = FALSE)
    
    if (any(lower.limits > 0)) {
      stop("Lower limits should be non-positive")
    }
    if (any(upper.limits < 0)) {
      stop("Upper limits should be non-negative")
    }
    lower.limits[lower.limits == -Inf] <- -control$big
    upper.limits[upper.limits == Inf] <- control$big
    if (length(lower.limits) < nvars) {
      if (length(lower.limits) == 1) {
        lower.limits <- rep(lower.limits, nvars)
      } else {
        stop("Require length 1 or nvars lower.limits")
      }
    } else {
      lower.limits <- lower.limits[seq(nvars)]
    }
    if (length(upper.limits) < nvars) {
      if (length(upper.limits) == 1) {
        upper.limits = rep(upper.limits, nvars)
      } else {
        stop("Require length 1 or nvars upper.limits")
      }
    } else {
      upper.limits <- upper.limits[seq(nvars)]
    }
    
    if (any(lower.limits == 0) || any(upper.limits == 0)) {
      ###Bounds of zero can mess with the lambda sequence and fdev;
      ###ie nothing happens and if fdev is not zero, the path can stop
      fdev <- multiview.control()$fdev
      if(fdev!= 0) {
        multiview.control(fdev = 0)
        on.exit(multiview.control(fdev = fdev))
      }
    }
    ### end check on limits
    ### end preparation of generic arguments
    # standardize x if necessary
    
    if (intercept) {
      xm <- col_means
      xm_list <- col_means_list
    } else {
      xm <- rep(0.0, times = nvars)
      xm_list<- lapply(p_x, numeric)
    }
    if (standardize) {
      xs <- col_sd
      xs_list <- col_sd_list
    } else {
      xs <- rep(1.0, times = nvars)
      xs_list <- lapply(p_x, rep, x = 1.0)
    }
    if (any(sapply(x_list, inherits, what = "sparseMatrix"))) {
      ## TODO: Something else needs to be done for SPARSE MATRIX
      nx_list <- x_list
      x <- do.call(cbind, x_list)
      attr(x, "xm") <- xm
      attr(x, "xs") <- xs
    } else {
      nx_list <- mapply(scale, x_list, xm_list, xs_list, SIMPLIFY = FALSE)
      x <- do.call(cbind, nx_list)
    }
    lower.limits <- lower.limits * xs
    upper.limits <- upper.limits * xs


    ## See whether this is a call to cox or not
    ## Better way is to check for a family object... see ?family
    if (cox_case) {
      ## fit <- multiview.cox.path(x_list = x_list, y = y, rho = rho, weights = weights, lambda = lambda,
      ##                           offset = offset, alpha = alpha, nlambda = nlambda,
      ##                           lambda.min.ratio = lambda.min.ratio,
      ##                           standardize = standardize, intercept = intercept,
      ##                           thresh = thresh, exclude = exclude, penalty.factor = penalty.factor,
      ##                           lower.limits = lower.limits, upper.limits = upper.limits, maxit = maxit,
      ##                           trace.it = 0, x = x)
        fit <- multiview.cox.path(x_list = nx_list, x = x, y = y, rho = rho, weights = weights,
                              lambda = lambda, nlambda = nlambda,
                              lambda.min.ratio = lambda.min.ratio,
                              alpha = alpha, offset = offset,
                              standardize = standardize, intercept = intercept,
                              thresh = thresh, maxit = maxit, penalty.factor = vp,
                              exclude = exclude,
                              lower.limits = lower.limits,
                              upper.limits = upper.limits,
                              trace.it = trace.it,
                              nvars = nvars, nobs = nobs, xm = xm, xs = xs, control = control, vp = vp, vnames = vnames,
                              is.offset = is.offset)

      fit$mvlambda <- fit$lambda  ## Iffy!      
    } else { ##

      # get null deviance and lambda max
      start_val <- get_start(x, y, weights, family, intercept, is.offset,
                             offset, exclude, vp, alpha)
      
      # work out lambda values
      nlam <- as.integer(nlambda)
      user_lambda <- FALSE   # did user provide their own lambda values?
      if (is.null(lambda)) {
        if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
        
        # compute lambda max: to add code here
        lambda_max <- start_val$lambda_max
        
        # compute lambda sequence
        ulam <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                        length.out = nlam))
      } else {  # user provided lambda values
        user_lambda <- TRUE
        if (any(lambda < 0)) stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- length(lambda)
      }
      
      if (gaussian_case) {
        ## Do a direct call
        ## Assume dense matrix first!
        ymean <- mean(y)
        ends  <- cumsum(p_x)
        starts  <- c(1, ends[-nviews] + 1)
        beta_indices <- mapply(seq.int, starts, ends, SIMPLIFY = FALSE)
        pairs <- apply(utils::combn(nviews, 2), 2, identity, simplify = FALSE)
        npairs <- length(pairs)
        ## The next line below is suboptimal for sparse matrices
        ##nx_list <- lapply(split(seq_len(nvars), rep(seq_along(p_x), p_x)), function(ind) x[, ind])
        rows <- lapply(pairs, make_row, x_list = nx_list, p_x = p_x, rho = rho )
        xt <- do.call(rbind, c(list(x), rows))
        yt <- c(y - ymean, rep(0, length(pairs) * nobs))
        ## Weights, offsets have to be handled in a special way because we augment the x matrix
        ## We have to ensure it is of the right length!
        if (!is.null(weights)) {
          weights <- c(weights, rep(weights, npairs))
        }
        if (is.offset) {
          offset <- c(offset, rep(offset, npairs))
        }
        fit <- glmnet::glmnet(x = xt, y = yt, family = "gaussian", weights = weights, lambda = ulam, 
                              offset = offset, alpha = alpha, nlambda = nlam,
                              lambda.min.ratio = lambda.min.ratio,
                              standardize = FALSE, intercept = FALSE, maxit = maxit,
                              thresh = thresh, exclude = exclude, penalty.factor = penalty.factor,
                              lower.limits = lower.limits, upper.limits = upper.limits,
                              trace.it = trace.it)
        ## Fix up intercept due to scaling of x 
        fit$beta <- fit$beta / xs
        fit$a0 <- ymean - colSums(as.matrix(fit$beta) * xm)
      } else {
        ## All other cases besided std gaussian
        ## arg_list  <- c(list(x_list = x_list, y = y, rho  = rho, family = family,
        ##                     exclude = exclude, mvlambda = mvlambda, x = x),
        ##                glmnet_args)
        ## fit <- do.call(multiview.path, arg_list)
        fit <- multiview.path(x_list = nx_list, x = x, y = y, rho = rho, weights = weights,
                              lambda = ulam, nlambda = nlam,
                              user_lambda = user_lambda, 
                              ## lambda.min.ratio = lambda.min.ratio, ## Not needed any more
                              alpha = alpha, offset = offset, family = family,
                              standardize = standardize, intercept = intercept,
                              thresh = thresh, maxit = maxit, penalty.factor = vp,
                              exclude = exclude,
                              lower.limits = lower.limits,
                              upper.limits = upper.limits,
                              trace.it = trace.it,
                              nvars = nvars, nobs = nobs, xm = xm, xs = xs, control = control, vp = vp, vnames = vnames,
                              start_val = start_val, is.offset = is.offset)

        }
      fit$mvlambda <- ulam ## TODO: CHECK
    }
  } else {

    ## rho <= 0 case, or early fusion.

    ## Prepare to call glmnet
    x <- do.call(cbind, x_list)
    col_sd <- apply(x, 2, sd)
    
    ## We dispatch to the old gaussian for speed (one less iteration)

    if (gaussian_case) our_family <- "gaussian" else our_family <- family
    if (cox_case) {
      fit <- glmnet::glmnet(x = x, y = y, family = our_family, weights = weights, offset = offset,
                            alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                            lambda = lambda, standardize = standardize, ## intercept = intercept, ## NO intercept for Cox
                            thresh = thresh, maxit = maxit, penalty.factor = penalty.factor,
                            exclude = exclude, lower.limits = lower.limits,
                            upper.limits = upper.limits, trace.it = trace.it)                           

    } else {
      fit <- glmnet::glmnet(x = x, y = y, family = our_family, weights = weights, offset = offset,
                            alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                            lambda = lambda, standardize = standardize, intercept = intercept,
                            thresh = thresh, maxit = maxit, penalty.factor = penalty.factor,
                            exclude = exclude, lower.limits = lower.limits,
                            upper.limits = upper.limits, trace.it = trace.it)
    }
    ## Ensure the mvlambda is updated to reflect value used in glmnet call
    fit$mvlambda <- fit$lambda
  }

  fit$call <- this.call
  fit$offset <- is.offset
  fit$p_x  <- p_x
  fit$colnames_list  <- colnames_list
  fit$rho <- rho
  fit$family <- family
  fit$xsd <- col_sd
  class(fit) <- c("multiview", class(fit))
  if (!is.character(family) && family$family == "binomial") {
    ## set class names
    fit$classnames <- names(table(y))
  }
  fit
}

#' Plot coefficients from a "multiview" object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' `"multiview"` object. The paths are colored by the data views, from which the features come.
#'
#' @param x A fitted `"multiview"` model.
#' @param col_palette A set of colors to use for indicating different views. If `NULL`,
#' the function will use the color palette "Set1" from the `RColorBrewer` package.
#' @param label If `TRUE`, label the curves with variable sequence.
#' numbers.
#' @param \dots Other graphical parameters to plot.
#' @return a `NULL` value as this function is really meant for its side-effect of generating a plot.
#' @examples
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = matrix(rnorm(100 * 10), 100, 10)
#' y = rnorm(100)
#' fit1 = multiview(list(x=x,z=z), y, rho = 0)
#' plot(fit1, label = TRUE)
#' 
#' # Binomial
#' by = sample(c(0,1), 100, replace = TRUE)
#' fit2 = multiview(list(x=x,z=z), by, family = binomial(), rho=0.5)
#' plot(fit2, label=FALSE)
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = multiview(list(x=x,z=z), py, family = poisson(), rho=0.5)
#' plot(fit3, label=TRUE)
#'
#' @method plot multiview
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats approx
#' @importFrom graphics axis legend matplot text
#' @export
plot.multiview  <- function(x,
                            col_palette = NULL, #rainbow(6, s = 0.5),
                            label = FALSE, ...) {
  object <- x
  beta_list <- object$beta
  lambda_list <- object$lambda
  rho <- object$rho
  index <- apply(abs(beta_list), 2, sum)
  n_coef <- sapply(seq_along(lambda_list), FUN = function(i) sum(beta_list[, i] != 0))
  p_x <- object$p_x  ## the p's for each view x
  if (is.null(col_palette)){
    col_palette <- brewer.pal(max(length(p_x), 3L), "Set1")[seq_along(p_x)] #rainbow(6, s = 0.5)
  }
  col_palette <- col_palette[seq_along(p_x)]
  color_list <- mapply(FUN = rep, col_palette, unlist(p_x))

  matplot(x = index, y = t(as.matrix(beta_list)),
          lty = 1, type = "l",
          xlab = 'L1 Norm', ylab = "Coefficients",
          col = unlist(color_list, use.names = FALSE),
          ...)
  atdf <- pretty(index)
  prettydf <- approx(x = index, y = n_coef, xout = atdf, rule = 2, method = "constant", f = 1)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  legend(x = 'topleft', legend = names(p_x), col = col_palette, lty = 1, bty = 'n')

  nr <- nrow(beta_list)
  beta <- abs(beta_list) > 0 # this is sparse
  which <- seq(nr)
  ones <- rep(1, ncol(beta_list))
  nz <- as.vector((beta_list %*% ones) > 0)
  which <- which[nz]
  #which=nonzeroCoef(beta)

  if (label){
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    xpos <- rep(xpos, nnz)
    ypos <- beta_list[, ncol(beta_list)]
    text(xpos, ypos, paste(which), cex = .5, pos = pos)
  }
}

#' Extract coefficients from a multiview object
#'
#' @examples
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = matrix(rnorm(100 * 10), 100, 10)
#' y = rnorm(100)
#' fit1 = multiview(list(x=x,z=z), y, rho = 0)
#' coef(fit1, s=0.1)
#' 
#' # Binomial
#' by = sample(c(0,1), 100, replace = TRUE)
#' fit2 = multiview(list(x=x,z=z), by, family = binomial(), rho=0.5)
#' coef(fit2, s=0.1)
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = multiview(list(x=x,z=z), py, family = poisson(), rho=0.5)
#' coef(fit3, s=0.1)
#' 
#' @method coef multiview
#' @inheritParams predict.multiview
#' @return a matrix of coefficients for specified lambda.
#' @export
coef.multiview <- function(object, s = NULL, ...) {
  #NextMethod(s = s, type = "coefficients", exact = exact, ...)
  predict(object, s = s, type = "coefficients", ...)
}

#' Extract an ordered list of standardized coefficients from a `multiview` or `cv.multiview` object
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
#' @param object Fitted `"multiview"` or `"cv.multiview"` object.
#'   coefficients are required.
#' @return data frame of consisting of view name, view column,
#'   coefficient and standardized coefficient ordered by rank of
#'   standardized coefficient.
#' @examples
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = matrix(rnorm(100 * 10), 100, 10)
#' y = rnorm(100)
#' fit1 = multiview(list(x=x,z=z), y, rho = 0)
#' coef_ordered(fit1, s=0.1)
#' 
#' # Binomial
#' by = sample(c(0,1), 100, replace = TRUE)
#' fit2 = multiview(list(x=x,z=z), by, family = binomial(), rho=0.5)
#' coef_ordered(fit2, s=0.1)
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = multiview(list(x=x,z=z), py, family = poisson(), rho=0.5)
#' coef_ordered(fit3, s=0.1)
#' @inheritParams predict.multiview
#' @export coef_ordered
coef_ordered <- function(object, ...) {
  UseMethod("coef_ordered")
}

## #' @export
## coef_ordered.default <- function(object, s = NULL, ...) {
##   stop(sprintf("coef_ordered method not defined for object of class: %s", class(object)))
## }

#' Extract an ordered list of standardized coefficients from a multiview object
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
#' @param object Fitted `"multiview"` object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#' coefficients are required.
#' @return data frame of consisting of view name, view column,
#'   coefficient and standardized coefficient ordered by rank of
#'   standardized coefficient.
#' @examples
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = matrix(rnorm(100 * 10), 100, 10)
#' y = rnorm(100)
#' fit1 = multiview(list(x=x,z=z), y, rho = 0)
#' coef_ordered(fit1, s=0.1)
#' 
#' # Binomial
#' by = sample(c(0,1), 100, replace = TRUE)
#' fit2 = multiview(list(x=x,z=z), by, family = binomial(), rho=0.5)
#' coef_ordered(fit2, s=0.1)
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = multiview(list(x=x,z=z), py, family = poisson(), rho=0.5)
#' coef_ordered(fit3, s=0.1)
#'
#' @method coef_ordered multiview
#' @inheritParams predict.multiview
#' @export
coef_ordered.multiview <- function(object, s = NULL, ...){
  if ((length(object$lambda) > 1) && (length(s) != 1)) {
    stop("s has to be specified as a single value")
  }
  beta <- as.numeric(coef(object, s = s, ...))[-1]
  which_nonzero <- which(beta != 0)
  p_x <- object$p_x
  view_names <- names(p_x)

  view <- c(rep(view_names, p_x))
  col_names  <- unlist(object$colnames_list)
  beta_std <- beta / object$xsd

  coefs  <- data.frame(view = view, view_col = col_names,
                       standardized_coef = beta_std, coef = beta)
  coefs_nonzero <- coefs[which_nonzero,]
  #coefs[rev(order(abs(beta_std))), ]
  coefs_nonzero[rev(order(abs(beta_std[which_nonzero]))), ]
}

process_exclude <- function(exclude, p_x) {
  if (!is.list(exclude) && (length(exclude) != length(p_x)) ||
        (!all(mapply(function(p, ind) all(1 <= ind & ind <= p), p_x, exclude)))) {
    stop("multiview: exclude list components do not match matrix dimensions in x_list")
  }
  exclude  <- lapply(exclude, unique)
  p_x  <- unlist(p_x)
  m  <- length(p_x)
  ends  <- cumsum(c(0L, p_x[-m]))
  exclude <- unlist(mapply(function(x, y) x + y, exclude, ends, SIMPLIFY = FALSE))

  if(length(exclude) > (sum(p_x) - 2L)) {
    stop("cannot retain 1 or less variables")
  }

  if(length(exclude) == 0) exclude <- NULL
  exclude
}

