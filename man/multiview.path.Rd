% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiview.path.R
\name{multiview.path}
\alias{multiview.path}
\title{Fit a GLM with elastic net regularization for a path of lambda values}
\usage{
multiview.path(
  x_list,
  y,
  rho = 0,
  weights = NULL,
  lambda,
  nlambda,
  user_lambda = FALSE,
  alpha = 1,
  offset = NULL,
  family = gaussian(),
  standardize = TRUE,
  intercept = TRUE,
  thresh = 1e-07,
  maxit = 1e+05,
  penalty.factor = rep(1, nvars),
  exclude = integer(0),
  lower.limits = -Inf,
  upper.limits = Inf,
  trace.it = 0,
  x,
  nvars,
  nobs,
  xm,
  xs,
  control,
  vp,
  vnames,
  start_val,
  is.offset
)
}
\arguments{
\item{x_list}{a list of \code{x} matrices with same number of rows
\code{nobs}}

\item{y}{the quantitative response with length equal to \code{nobs}, the
(same) number of rows in each \code{x} matrix}

\item{rho}{the weight on the agreement penalty, default 0. \code{rho=0}
is a form of early fusion, and \code{rho=1} is a form of late fusion.
We recommend trying a few values of \code{rho} including 0, 0.1, 0.25,
0.5, and 1 first; sometimes \code{rho} larger than 1 can also be
helpful.}

\item{weights}{observation weights. Can be total counts if
responses are proportion matrices. Default is 1 for each
observation}

\item{lambda}{A user supplied \code{lambda} sequence, default
\code{NULL}. Typical usage is to have the program compute its own
\code{lambda} sequence. This sequence, in general, is different from
that used in the \code{\link[glmnet:glmnet]{glmnet::glmnet()}} call (named \code{lambda})
Supplying a value of \code{lambda} overrides this. WARNING: use with
care. Avoid supplying a single value for \code{lambda} (for
predictions after CV use \code{\link[stats:predict]{stats::predict()}} instead.  Supply
instead a decreasing sequence of \code{lambda} values as \code{multiview}
relies on its warms starts for speed, and its often faster to fit
a whole path than compute a single fit.}

\item{nlambda}{The number of \code{lambda} values - default is 100.}

\item{user_lambda}{a flag indicating if user supplied the lambda sequence}

\item{alpha}{The elasticnet mixing parameter, with
\eqn{0\le\alpha\le 1}.  The penalty is defined as
\deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.}
\code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge
penalty.}

\item{offset}{A vector of length \code{nobs} that is included in
the linear predictor (a \code{nobs x nc} matrix for the
\code{"multinomial"} family).  Useful for the \code{"poisson"}
family (e.g. log of exposure time), or for refining a model by
starting at a current fit. Default is \code{NULL}. If supplied,
then values must also be supplied to the \code{predict} function.}

\item{family}{A description of the error distribution and link
function to be used in the model. This is the result of a call to
a family function. Default is \link[stats:family]{stats::gaussian}. (See
\link[stats:family]{stats::family} for details on family functions.)}

\item{standardize}{Logical flag for x variable standardization,
prior to fitting the model sequence. The coefficients are always
returned on the original scale. Default is
\code{standardize=TRUE}.  If variables are in the same units
already, you might not wish to standardize. See details below for
y standardization with \code{family="gaussian"}.}

\item{intercept}{Should intercept(s) be fitted (default \code{TRUE})}

\item{thresh}{Convergence threshold for coordinate descent. Each
inner coordinate-descent loop continues until the maximum change
in the objective after any coefficient update is less than
\code{thresh} times the null deviance. Defaults value is
\code{1E-7}.}

\item{maxit}{Maximum number of passes over the data for all lambda
values; default is 10^5.}

\item{penalty.factor}{Separate penalty factors can be applied to
each coefficient. This is a number that multiplies \code{lambda}
to allow differential shrinkage. Can be 0 for some variables,
which implies no shrinkage, and that variable is always included
in the model. Default is 1 for all variables (and implicitly
infinity for variables listed in \code{exclude}). Note: the
penalty factors are internally rescaled to sum to nvars, and the
lambda sequence will reflect this change.}

\item{exclude}{Indices of variables to be excluded from the
model. Default is none. Equivalent to an infinite penalty factor
for the variables excluded (next item).  Users can supply instead
an \code{exclude} function that generates the list of indices.  This
function is most generally defined as \verb{function(x_list, y, ...)},
and is called inside \code{multiview} to generate the indices for
excluded variables.  The \code{...} argument is required, the others
are optional.  This is useful for filtering wide data, and works
correctly with \code{cv.multiview}. See the vignette 'Introduction'
for examples.}

\item{lower.limits}{Vector of lower limits for each coefficient;
default \code{-Inf}. Each of these must be non-positive. Can be
presented as a single value (which will then be replicated), else
a vector of length \code{nvars}}

\item{upper.limits}{Vector of upper limits for each coefficient;
default \code{Inf}. See \code{lower.limits}}

\item{trace.it}{If \code{trace.it=1}, then a progress bar is
displayed; useful for big models that take a long time to fit.}

\item{x}{the \code{cbind}ed matrices in \code{x_list}}

\item{nvars}{the number of variables (total)}

\item{nobs}{the number of observations}

\item{xm}{the column means vector (could be zeros if \code{standardize = FALSE})}

\item{xs}{the column std dev vector (could be 1s if \code{standardize = FALSE})}

\item{control}{the multiview control object}

\item{vp}{the variable penalities (processed)}

\item{vnames}{the variable names}

\item{start_val}{the result of first call to \code{get_start}}

\item{is.offset}{a flag indicating if offset is supplied or not}
}
\value{
An object with class \code{"multiview"} \code{"glmnetfit"} and \code{"glmnet"}
\item{a0}{Intercept sequence of length \code{length(lambda)}.}
\item{beta}{A \verb{nvars x length(lambda)} matrix of coefficients, stored in
sparse matrix format.}
\item{df}{The number of nonzero coefficients for each value of lambda.}
\item{dim}{Dimension of coefficient matrix.}
\item{lambda}{The actual sequence of lambda values used. When alpha=0, the
largest lambda reported does not quite give the zero coefficients reported
(lambda=inf would in principle). Instead, the largest lambda for alpha=0.001
is used, and the sequence of lambda values is derived from this.}
\item{lambda}{The sequence of lambda values}
\item{mvlambda}{The corresponding sequence of multiview lambda values}
\item{dev.ratio}{The fraction of (null) deviance explained. The deviance
calculations incorporate weights if present in the model. The deviance is
defined to be 2*(loglike_sat - loglike), where loglike_sat is the log-likelihood
for the saturated model (a model with a free parameter per observation).
Hence dev.ratio=1-dev/nulldev.}
\item{nulldev}{Null deviance (per observation). This is defined to be
2*(loglike_sat -loglike(Null)). The null model refers to the intercept model.}
\item{npasses}{Total passes over the data summed over all lambda values.}
\item{jerr}{Error flag, for warnings and errors (largely for internal
debugging).}
\item{offset}{A logical variable indicating whether an offset was included
in the model.}
\item{call}{The call that produced this object.}
\item{family}{Family used for the model.}
\item{nobs}{Number of observations.}
}
\description{
Fit a generalized linear model via penalized maximum likelihood for a path of
lambda values. Can deal with any GLM family.
}
\details{
\code{multiview.path} solves the elastic net problem for a path of lambda values.
It generalizes \code{multiview::multiview} in that it works for any GLM family.

Sometimes the sequence is truncated before \code{nlam} values of lambda
have been used. This happens when \code{multiview.path} detects that the decrease
in deviance is marginal (i.e. we are near a saturated fit).
}
