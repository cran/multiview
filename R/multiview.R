#' Perform cooperative learning using the direct algorithm for
#' two or more views.
#'
#' `multiview` uses [glmnet::glmnet()] to do most of its work and
#' therefore takes many of the same parameters, but an intercept is
#' always included, standardization is always done and several other
#' parameters do not apply. Therefore they are always overridden and
#' warnings issued.
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
#' @param x_list a list of `x` matrices with same number of rows
#'   `nobs`
#' @param y the quantitative response with length equal to `nobs`, the
#'   (same) number of rows in each `x` matrix
#' @param rho the weight on the agreement penalty, default 0. `rho=0`
#'   is a form of early fusion, and `rho=1` is a form of late fusion.
#'   We recommend trying a few values of `rho` including 0, 0.1, 0.25, 0.5, and 1 first; 
#'   sometimes `rho` larger than 1 can also be helpful.
#' @param family A description of the error distribution and link
#'   function to be used in the model. This is the result of a call to
#'   a family function. Default is [stats::gaussian]. (See
#'   [stats::family] for details on family functions.)
#' @param exclude Indices of variables to be excluded from the
#'   model. Default is none. Equivalent to an infinite penalty factor
#'   for the variables excluded (next item).  Users can supply instead
#'   an `exclude` function that generates the list of indices.  This
#'   function is most generally defined as `function(x_list, y, ...)`,
#'   and is called inside `multiview` to generate the indices for
#'   excluded variables.  The `...` argument is required, the others
#'   are optional.  This is useful for filtering wide data, and works
#'   correctly with `cv.glmnet`.
#' @param mvlambda A user supplied `lambda` sequence, default
#'   `NULL`. Typical usage is to have the program compute its own
#'   `mvlambda` sequence. This sequence, in general, is different from
#'   that used in the [glmnet::glmnet()] call (named `lambda`)
#'   Supplying a value of `mvlambda` overrides this. WARNING: use with
#'   care. Avoid supplying a single value for `mvlambda` (for
#'   predictions after CV use [stats::predict()] instead.  Supply
#'   instead a decreasing sequence of `mvlambda` values as `multiview`
#'   relies on its warms starts for speed, and its often faster to fit
#'   a whole path than compute a single fit.
#' @param ... further arguments to glmnet, some of which may be
#'   overridden as noted above
#' @return An object with S3 class `"multiview","*" `, where `"*"` is
#'   `"elnet"`, `"lognet"`, `"multnet"`, `"fishnet"` (poisson),
#'   `"coxnet"` or `"mrelnet"` for the various types of models.
#'   \item{call}{the call that produced this object}
#'   \item{a0}{Intercept sequence of length `length(lambda)`}
#'   \item{beta}{For `"elnet"`, `"lognet"`, `"fishnet"` and `"coxnet"`
#'   models, a `nvars x length(lambda)` matrix of coefficients, stored
#'   in sparse column format (`"CsparseMatrix"`). For `"multnet"` and
#'   `"mgaussian"`, a list of `nc` such matrices, one for each class.}
#'   \item{lambda}{The actual sequence of [glmnet::glmnet()] `lambda` values used. When
#'   `alpha=0`, the largest lambda reported does not quite give the
#'   zero coefficients reported (`lambda=inf` would in principle).
#'   Instead, the largest `lambda` for `alpha=0.001` is used, and the
#'   sequence of `lambda` values is derived from this.}
#' \item{mvlambda}{The corresponding sequence of multiview lambda values}
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
multiview <- function(x_list, y, rho = 0, family = gaussian(), 
                      exclude = NULL, mvlambda = NULL, ...) {
  this.call <- match.call()

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

  ## Prepare to reuse glmnetFlex code
  x <- do.call(cbind, x_list)
  ## We need the std devs for other purposes, so we compute it
  xsd <- apply(x, 2, sd)

  ## Handle exclude like in glmnet
  if(is.function(exclude)) exclude <- process.exclude(exclude(x_list = x_list, y = y, ...), p_x)

  ## Ensure that some arguments are ignored or kept at default values for glmnet
  glmnet_args  <- list(...)
  ## Use of the args below will emit a warning
  forced_args  <- list(intercept = TRUE, standardize = TRUE, standardize.response = FALSE, relax = FALSE, alpha = 1.0)
  specified_args <- intersect(names(glmnet_args), names(forced_args))
  if (length(specified_args) > 0) {
    for (arg in specified_args) {
      if (glmnet_args[[arg]] != forced_args[[arg]]) {
        glmnet_args[[arg]] <- forced_args[[arg]]
        warning(sprintf("multiview: overriding specified %s to %s", arg, forced_args[[arg]]))
      }
    }
  }

  ## Glmnet lambda is never specifiable.
  if (!is.null(glmnet_args$lambda)) {
    warning("multiview: glmnet lambda is computed internally, so ignored")
    glmnet_args$lambda <- NULL
  }

  ## Internally, of course, we standardize, but our call to glmnet
  ## needs to have standardize = FALSE. Yes, this is bloody confusing
  ## and we can go mad if we muck this up!
  glmnet_args$standardize <- FALSE

  ## direct_call <- (family$family == "gaussian" && family$link == "identity")
  ## if (direct_call) {
  ##   arg_list  <- c(list(x = x, y = y, rho  = rho, family = family),
  ##                  glmnet_args)
  ##   fit <- do.call(glmnet::glmnet, arg_list)
  ##   ## Ensure the mvlambda is updated to reflect value used in glmnet call
  ##   fit$mvlambda <- fit$lambda
  ## }

  gaussian_case <- (family$family == "gaussian" && family$link == "identity")
  
  if (rho > 0 && nviews > 1L) {

    ## See whether this is a call to cox or not
    ## Better way is to check for a family object... see ?family
    if (is.character(family)) {
      if (family != "cox") {
        stop("multiview: family must be a family function or the string 'cox'")
      }
      ## fit <- multiview.cox.path(x_list, y, weights, offset, alpha, nlambda, lambda.min.ratio,
      ##                           lambda, standardize, thresh, exclude, penalty.factor,
      ##                           lower.limits, upper.limits, maxit, trace.it, ...)
      stop("multiview: Penalized Cox model not implemented yet!")
    } else { ##
      if (gaussian_case) {
        ## Do a direct call
        nx_list <- lapply(x_list, scale, center = TRUE, scale = FALSE)
        xm <- unlist(lapply(nx_list, function(mat) attr(mat, "scaled:center")))
        nx <- do.call(cbind, nx_list)
        ymean <- mean(y)
        ends  <- cumsum(p_x)
        starts  <- c(1, ends[-nviews] + 1)
        beta_indices <- mapply(seq.int, starts, ends, SIMPLIFY = FALSE)
        pairs <- apply(utils::combn(nviews, 2), 2, identity, simplify = FALSE)
        npairs <- length(pairs)
        rows <- lapply(pairs, make_row, x_list = nx_list, p_x = p_x, rho = rho )
        xt <- do.call(rbind, c(list(nx), rows))
        ## yt <- c(y, rep(0, length(pairs) * nobs))
        yt <- c(y - ymean, rep(0, length(pairs) * nobs))
        glmnet_args$standardize <- FALSE
        glmnet_args$intercept <- FALSE
        ## Weights have to be handled in a special way because we augment the x matrix
        ## We have to ensure it is of the right length!
        if (!is.null(glmnet_args$weights)) {
          glmnet_args$weights <- c(glmnet_args$weights, rep(glmnet_args$weights, npairs))
        }
        arg_list <- c(list(x = xt, y = yt, exclude = exclude, family = "gaussian"),
                      glmnet_args)
        fit <- do.call(glmnet::glmnet, arg_list)
        ## Fix up intercept due to scaling of x 
        ## fit$a0 <- fit$a0 - colSums(as.matrix(fit$beta) * xm)
        ## fit$a0 <- fit$a0 * (npairs + 1)
        fit$a0 <- ymean - colSums(as.matrix(fit$beta) * xm)
      } else {
        arg_list  <- c(list(x_list = x_list, y = y, rho  = rho, family = family,
                            exclude = exclude, mvlambda = mvlambda, x = x),
                       glmnet_args)
        fit <- do.call(multiview.path, arg_list)
      }
      fit$mvlambda <- fit$lambda  ## Iffy!      
    }
  } else {
    ## rho <= 0 case, or early fusion.
    ## We dispatch to the old gaussian for speed (one less iteration)
    our_family <- family
    if (family$family == "gaussian" && family$link == "identity") our_family <- "gaussian"
    glmnet_args$lambda <- mvlambda
    arg_list  <- c(list(x = x, y = y, rho  = rho, family = our_family,
                        exclude = exclude),
                   glmnet_args)
    fit <- do.call(glmnet::glmnet, arg_list)
    ## Ensure the mvlambda is updated to reflect value used in glmnet call
    fit$mvlambda <- fit$lambda
  }

  fit$call <- this.call
  fit$p_x  <- p_x
  fit$colnames_list  <- colnames_list
  fit$rho <- rho
  fit$family <- family
  fit$xsd <- xsd
  class(fit) <- c("multiview", class(fit))
  if (family$family == "binomial") {
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

process.exclude <- function(exclude, p_x) {
  if ((length(exclude) != length(p_x)) ||
        (!all(mapply(function(p, ind) all(1 <= ind & ind <= p), p_x, exclude)))) {
    stop("multiview: exclude should be list of column indices (possibly NULL) for each matrix in x_list")
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

