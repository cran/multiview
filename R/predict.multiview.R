#' Get predictions from a `multiview` fit object
#'
#' Gives fitted values, linear predictors, coefficients and number of non-zero
#' coefficients from a fitted `multiview` object.
#'
#' @param object Fitted "multiview" object.
#' @param newx list of new matrices for  `x` at which predictions are to be
#' made. Must be a list of matrices. This argument is not used for `type =
#' c("coefficients","nonzero")`.
#' @param s Value(s) of the penalty parameter lambda at which predictions are
#' required. Default is the entire sequence used to create the model.
#' @param type Type of prediction required. Type "link" gives the linear
#' predictors (eta scale); Type "response" gives the fitted values (mu scale).
#' Type "coefficients" computes the coefficients at the requested values for s.
#' Type "nonzero" returns a list of the indices of the nonzero coefficients for
#' each value of s. Type "class" returns class labels for binomial family only. 
#' @param exact This argument is relevant only when predictions are made at values
#' of `s` (`lambda`) \emph{different} from those used in the fitting of the
#' original model. If `exact=FALSE` (default), then the predict function
#' uses linear interpolation to make predictions for values of `s` (lambda)
#' that do not coincide with those used in the fitting algorithm. While this is
#' often a good approximation, it can sometimes be a bit coarse. With
#' `exact=TRUE`, these different values of `s` are merged (and sorted)
#' with `object$lambda`, and the model is refit before predictions are made.
#' In this case, it is required to supply the original data x= and y= as additional
#' named arguments to predict() or coef(). The workhorse `predict.multiview()`
#' needs to update the model, and so needs the data used to create it. The same
#' is true of weights, offset, penalty.factor, lower.limits, upper.limits if
#' these were used in the original call. Failure to do so will result in an error.
#' @param newoffset If an offset is used in the fit, then one must be supplied for
#' making predictions (except for type="coefficients" or type="nonzero").
#' @param ... This is the mechanism for passing arguments like `x=` when
#' `exact=TRUE`; see `exact` argument.
#'
#' @return The object returned depends on type.
#' 
#' 
#' @examples
#' # Gaussian
#' x = matrix(rnorm(100 * 20), 100, 20)
#' z = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = multiview(list(x=x,z=z), y, rho = 0)
#' predict(fit1, newx = list(x[1:10, ],z[1:10, ]), s = c(0.01, 0.005))
#' 
#' # Binomial
#' by = sample(c(0,1), 100, replace = TRUE)
#' fit2 = multiview(list(x=x,z=z), by, family = binomial(), rho=0.5)
#' predict(fit2, newx = list(x[1:10, ],z[1:10, ]), s = c(0.01, 0.005), type = "response")
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = multiview(list(x=x,z=z), py, family = poisson(), rho=0.5)
#' predict(fit3, newx = list(x[1:10, ],z[1:10, ]), s = c(0.01, 0.005), type = "response")
#'
#' @method predict multiview
#' @export
predict.multiview <- function(object, newx, s = NULL,
                              type = c("link", "response", "coefficients",  "class", "nonzero"),
                              exact = FALSE, newoffset, ...) {
  ## type <- match.arg(type)
  ##   if (!missing(newx)) newx <- do.call(cbind, newx)
  ##   nfit <- NextMethod("predict")
  ##   if (type == "response") {
  ##     object$family$linkinv(nfit)
  ##   } else {
  ##       nfit
  ##   }

  type <- match.arg(type)


  if (type %in% c("coefficients", "nonzero")) {
    return(NextMethod(s = s, type = type, ...))
  }

  newxt <- do.call(cbind, newx)

  type_class <- type == "class"
  if (type_class) {
    if (object$family$family != "binomial") {
      stop("predict.multiview: type = 'class' is only for binomial family")
    }
    type <- "link"  # modify for type == 'class' only.
  }

  if (missing(newoffset)) {
    nfit <- NextMethod(newx = newxt, s = s, type = type, exact = exact, ...)
  } else {
    nfit <- NextMethod(newx = newxt, s = s, type = type, exact = exact, newoffset = newoffset, ...)
  }

  if (type_class) {
    cnum <- ifelse(nfit > 0, 2, 1)
    clet <- object$classnames[cnum]
    if (is.matrix(cnum)) clet <- array(clet, dim(cnum), dimnames(cnum))
    clet
  } else {
    nfit
  }
}

