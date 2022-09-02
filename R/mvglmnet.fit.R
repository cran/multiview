# Copy of glmnet.fit for use in get_start.R See lines 57 onwards where if some penalty factors
# are zero, mu has to be recomputed
mvglmnet.fit <- function(x, y, weights, lambda, alpha = 1.0,
                       offset = rep(0, nobs), family = gaussian(),
                       intercept = TRUE, thresh = 1e-10, maxit = 100000,
                       penalty.factor = rep(1.0, nvars), exclude = c(), lower.limits = -Inf,
                       upper.limits = Inf, warm = NULL, from.glmnet.path = FALSE,
                       save.fit = FALSE, trace.it = 0) {
    this.call <- match.call()
    control <- multiview.control()

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
    if (!from.glmnet.path) {
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
    lower.limits[lower.limits == -Inf] = -control$big
    upper.limits[upper.limits == Inf] = control$big
    if (length(lower.limits) < nvars)
        lower.limits = rep(lower.limits, nvars) else
            lower.limits = lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars)
        upper.limits = rep(upper.limits, nvars) else
            upper.limits = upper.limits[seq(nvars)]
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
    obj_val_old <- obj_function(y, mu, weights, family, lambda, alpha, coefold, vp)
    if (trace.it == 2) {
        cat("Warm Start Objective:", obj_val_old, fill = TRUE)
    }
    conv <- FALSE      # converged?

    # IRLS loop
    for (iter in 1L:control$mxitnr) {
        # some checks for NAs/zeros
        varmu <- variance(mu)
        if (anyNA(varmu)) stop("NAs in V(mu)")
        if (any(varmu == 0)) stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")

        # compute working response and weights
        z <- (eta - offset) + (y - mu)/mu.eta.val
        w <- (weights * mu.eta.val^2)/variance(mu)

        # have to update the weighted residual in our fit object
        # (in theory g and iy should be updated too, but we actually recompute g
        # and iy anyway in wls.f)
        if (!is.null(fit)) {
            fit$warm_fit$r <- w * (z - eta + offset)
        }

        # do WLS with warmstart from previous iteration
        fit <- elnet.fit(x, z, w, lambda, alpha, intercept,
                         thresh = thresh, maxit = maxit, penalty.factor = vp,
                         exclude = exclude, lower.limits = lower.limits,
                         upper.limits = upper.limits, warm = fit,
                         from.glmnet.fit = TRUE, save.fit = TRUE)
      if (fit$jerr != 0) return(list(jerr = fit$jerr))

        # update coefficients, eta, mu and obj_val
        start <- fit$warm_fit$a
        start_int <- fit$warm_fit$aint
        eta <- get_eta(x, start, start_int)
        mu <- linkinv(eta <- eta + offset)
        obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
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
                obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
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
            obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
            if (trace.it == 2) cat("Iteration", iter, " Halved step 2, Objective:", obj_val, fill = TRUE)
        }
        # extra halving step if objective function value actually increased
        if (obj_val > obj_val_old + 1e-7) {
            ii <- 1
            while (obj_val > obj_val_old + 1e-7) {
                if (ii > control$mxitnr)
                    stop("inner loop 3; cannot correct step size", call. = FALSE)
                ii <- ii + 1
                start <- (start + coefold)/2
                start_int <- (start_int + intold)/2
                eta <- get_eta(x, start, start_int)
                mu <- linkinv(eta <- eta + offset)
                obj_val <- obj_function(y, mu, weights, family, lambda, alpha, start, vp)
                if (trace.it == 2) cat("Iteration", iter, " Halved step 3, Objective:",
                                       obj_val, fill = TRUE)
            }
            halved <- TRUE
        }

        # if we did any halving, we have to update the coefficients, intercept
        # and weighted residual in the warm_fit object
        if (halved) {
            fit$warm_fit$a <- start
            fit$warm_fit$aint <- start_int
            fit$warm_fit$r <- w * (z - eta)
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

    # add new key-value pairs to list
    fit$family <- family
    fit$converged <- conv
    fit$boundary <- boundary
    fit$obj_function <- obj_val

    class(fit) <- c("glmnetfit", "glmnet")
    fit
}
