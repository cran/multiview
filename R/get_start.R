#' Get null deviance, starting mu and lambda max
#'
#' Return the null deviance, starting mu and lambda max values for
#' initialization. For internal use only.
#'
#' This function is called by \code{glmnet.path} for null deviance, starting mu
#' and lambda max values. It is also called by \code{glmnet.fit} when used
#' without warmstart, but they only use the null deviance and starting mu values.
#'
#' When \code{x} is not sparse, it is expected to already by centered and scaled.
#' When \code{x} is sparse, the function will get its attributes \code{xm} and
#' \code{xs} for its centering and scaling factors.
#'
#' Note that whether \code{x} is centered & scaled or not, the values of \code{mu}
#' and \code{nulldev} don't change. However, the value of \code{lambda_max} does
#' change, and we need \code{xm} and \code{xs} to get the correct value.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed to be standardized.
#' @param y Quantitative response variable.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
#' (See \code{\link[stats:family]{family}} for details on family functions.)
#' @param intercept Does the model we are fitting have an intercept term or not?
#' @param is.offset Is the model being fit with an offset or not?
#' @param offset Offset for the model. If \code{is.offset=FALSE}, this should be
#' a zero vector of the same length as \code{y}.
#' @param exclude Indices of variables to be excluded from the model.
#' @param vp Separate penalty factors can be applied to each coefficient.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @importFrom stats glm glm.fit
#' 
get_start <- function(x, y, weights, family, intercept, is.offset, offset,
                      exclude, vp, alpha) {
    nobs <- nrow(x); nvars <- ncol(x)

    # compute mu and null deviance
    # family = binomial() gives us warnings due to non-integer weights
    # to avoid, suppress warnings
    if (intercept) {
        if (is.offset) {
            suppressWarnings(tempfit <- glm(y ~ 1, family = family,
                                            weights = weights, offset = offset))
            mu <- tempfit$fitted.values
        } else {
            mu <- rep(weighted.mean(y,weights), times = nobs)
        }
    } else {
        mu <- family$linkinv(offset)
    }
    ## nulldev <- dev_function(y, mu, weights, family)
    nulldev <- dev_function(y, mu, weights, family)

    # if some penalty factors are zero, we have to recompute mu
    vp_zero <- setdiff(which(vp == 0), exclude)
    if (length(vp_zero) > 0) {
      tempx <- x[, vp_zero, drop = FALSE]
      if (inherits(tempx, "sparseMatrix")) {
        ## tempfit <- glmnet.fit(tempx, y, intercept = intercept, family = family,
        ##                       weights = weights/sum(weights), offset = offset, lambda = 0)
        tempfit <- mvglmnet.fit(tempx, y, intercept = intercept, family = family,
                                weights = weights/sum(weights), offset = offset, lambda = 0)
        mu <- predict(tempfit, newx=tempx, newoffset=offset, type = "response")
      } else {
        if (intercept) {
          tempx <- cbind(1,tempx)
        }
        tempfit <- glm.fit(tempx, y, family = family, weights = weights, offset = offset)
        mu <- tempfit$fitted.values
      }
    }
    # compute lambda max
    ju <- rep(1, nvars)
    ju[exclude] <- 0 # we have already included constant variables in exclude
    r <- y - mu
    eta <- family$linkfun(mu)
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)
    weights <- weights / sum(weights)
    rv <- r / v * m.e * weights
    if (inherits(x, "sparseMatrix")) {
        xm <- attr(x, "xm")
        xs <- attr(x, "xs")
        g <- abs((drop(t(rv) %*% x) - sum(rv) * xm) / xs)
    } else {
        g <- abs(drop(t(rv) %*% x))
    }
    g <- g * ju / ifelse(vp > 0, vp, 1)
    lambda_max <- max(g) / max(alpha, 1e-3)

    list(nulldev = nulldev, mu = mu, lambda_max = lambda_max)
}

#' Elastic net objective function value
#'
#' Returns the elastic net objective function value.
#'
#' @param y Quantitative response variable.
#' @param mu Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
#' @param lambda A single value for the \code{lambda} hyperparameter.
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param coefficients The model's coefficients (excluding intercept).
#' @param vp Penalty factors for each of the coefficients.
#' @param view_components a list of lists containing indices of coefficients and associated covariate (view) pairs
#' @param rho the fusion parameter
obj_function <- function(y, mu, weights, family,
                         lambda, alpha, coefficients, vp, view_components, rho) {

  coop_terms <- lapply(view_components, function(l) {
    sum((l$x[[1L]] %*% coefficients[ l$index[[1L]] ] -
          l$x[[2L]] %*% coefficients[ l$index[[2L]] ])^2)
  })

  dev_function(y, mu, weights, family) / 2  +
    #lambda * pen_function(coefficients, alpha, vp) +
    # inlining the definition to avoid a function call
    lambda * sum(vp * (alpha * abs(coefficients) + (1-alpha)/2 * coefficients^2)) +
    0.5 * rho *  Reduce(f = '+', x = coop_terms)
}

#' Elastic net penalty value
#'
#' Returns the elastic net penalty value without the \code{lambda} factor.
#'
#' The penalty is defined as
#' \deqn{(1-\alpha)/2 \sum vp_j \beta_j^2 + \alpha \sum vp_j |\beta|.}
#' Note the omission of the multiplicative \code{lambda} factor.
#'
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}.
#' @param coefficients The model's coefficients (excluding intercept).
#' @param vp Penalty factors for each of the coefficients.
pen_function <- function(coefficients, alpha = 1.0, vp = 1.0) {
    sum(vp * (alpha * abs(coefficients) + (1-alpha)/2 * coefficients^2))
}

#' Elastic net deviance value
#'
#' Returns the elastic net deviance value.
#'
#' @param y Quantitative response variable.
#' @param mu Model's predictions for \code{y}.
#' @param weights Observation weights.
#' @param family A description of the error distribution and link function to be
#' used in the model. This is the result of a call to a family function.
dev_function <- function(y, mu, weights, family) {
    sum(family$dev.resids(y, mu, weights))
}

jerr.glmnetfit <- function (n, maxit, k = NULL) {
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

#' Helper function to get etas (linear predictions)
#'
#' Given x, coefficients and intercept, return linear predictions. Wrapper that
#' works with both regular and sparse x. Only works for single set of coefficients
#' and intercept.
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is an
#' observation vector. If it is a sparse matrix, it is assumed to be unstandardized.
#' It should have attributes \code{xm} and \code{xs}, where \code{xm(j)} and
#' \code{xs(j)} are the centering and scaling factors for variable j respsectively.
#' If it is not a sparse matrix, it is assumed to be standardized.
#' @param beta Feature coefficients.
#' @param a0 Intercept.
get_eta <- function(x, beta, a0) {
    if (inherits(x, "sparseMatrix")) {
        beta <- drop(beta)/attr(x, "xs")
        as.numeric(x %*% beta - sum(beta * attr(x, "xm") ) + a0)
    } else {
        drop(x %*% beta + a0)
    }
}

#' Helper function to compute weighted mean and standard deviation
#'
#' Helper function to compute weighted mean and standard deviation.
#' Deals gracefully whether x is sparse matrix or not.
#'
#' @param x Observation matrix.
#' @param weights Optional weight vector.
#'
#' @return A list with components.
#' \item{mean}{vector of weighted means of columns of x}
#' \item{sd}{vector of weighted standard deviations of columns of x}
weighted_mean_sd <- function(x, weights=rep(1,nrow(x))){
    weights <- weights/sum(weights)
    xm <- drop(t(weights)%*%x)
    xv <- drop(t(weights)%*%scale(x,xm,FALSE)^2)
    xv[xv < 10*.Machine$double.eps] <- 0
    list(mean = xm, sd = sqrt(xv))
}

