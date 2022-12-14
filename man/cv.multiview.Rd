% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.multiview.R
\name{cv.multiview}
\alias{cv.multiview}
\title{Perform k-fold cross-validation for cooperative learning}
\usage{
cv.multiview(
  x_list,
  y,
  family = gaussian(),
  rho = 0,
  weights = NULL,
  offset = NULL,
  lambda = NULL,
  type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
  nfolds = 10,
  foldid = NULL,
  alignment = c("lambda", "fraction"),
  grouped = TRUE,
  keep = FALSE,
  trace.it = 0,
  ...
)
}
\arguments{
\item{x_list}{a list of \code{x} matrices with same number of rows
\code{nobs}}

\item{y}{the quantitative response with length equal to \code{nobs}, the
(same) number of rows in each \code{x} matrix}

\item{family}{A description of the error distribution and link
function to be used in the model. This is the result of a call to
a family function. Default is \link[stats:family]{stats::gaussian}. (See
\link[stats:family]{stats::family} for details on family functions.)}

\item{rho}{the weight on the agreement penalty, default 0. \code{rho=0}
is a form of early fusion, and \code{rho=1} is a form of late fusion.
We recommend trying a few values of \code{rho} including 0, 0.1, 0.25,
0.5, and 1 first; sometimes \code{rho} larger than 1 can also be
helpful.}

\item{weights}{Observation weights; defaults to 1 per observation}

\item{offset}{Offset vector (matrix) as in \code{multiview}}

\item{lambda}{A user supplied \code{lambda} sequence, default
\code{NULL}. Typical usage is to have the program compute its own
\code{lambda} sequence. This sequence, in general, is different from
that used in the \code{\link[glmnet:glmnet]{glmnet::glmnet()}} call (named \code{lambda}). Note
that this is done for the full model (master sequence), and
separately for each fold.  The fits are then aligned using the
glmnet lambda sequence associated with the master sequence (see
the \code{alignment} argument for additional details). Adapting
\code{lambda} for each fold leads to better convergence. When
\code{lambda} is supplied, the same sequence is used everywhere, but
in some GLMs can lead to convergence issues.}

\item{type.measure}{loss to use for cross-validation. Currently
five options, not all available for all models. The default is
\code{type.measure="deviance"}, which uses squared-error for gaussian
models (a.k.a \code{type.measure="mse"} there), deviance for logistic
and poisson regression, and partial-likelihood for the Cox model.
\code{type.measure="class"} applies to binomial and multinomial
logistic regression only, and gives misclassification error.
\code{type.measure="auc"} is for two-class logistic regression only,
and gives area under the ROC curve. \code{type.measure="mse"} or
\code{type.measure="mae"} (mean absolute error) can be used by all
models except the \code{"cox"}; they measure the deviation from the
fitted mean to the response.  \code{type.measure="C"} is Harrel's
concordance measure, only available for \code{cox} models.}

\item{nfolds}{number of folds - default is 10. Although \code{nfolds}
can be as large as the sample size (leave-one-out CV), it is not
recommended for large datasets. Smallest value allowable is
\code{nfolds=3}}

\item{foldid}{an optional vector of values between 1 and \code{nfold}
identifying what fold each observation is in. If supplied,
\code{nfold} can be missing.}

\item{alignment}{This is an experimental argument, designed to fix
the problems users were having with CV, with possible values
\code{"lambda"} (the default) else \code{"fraction"}. With \code{"lambda"} the
\code{lambda} values from the master fit (on all the data) are used to
line up the predictions from each of the folds. In some cases
this can give strange values, since the effective \code{lambda} values
in each fold could be quite different. With \code{"fraction"} we line
up the predictions in each fold according to the fraction of
progress along the regularization. If in the call a \code{lambda}
argument is also provided, \code{alignment="fraction"} is ignored
(with a warning).}

\item{grouped}{This is an experimental argument, with default
\code{TRUE}, and can be ignored by most users. For all models except
the \code{"cox"}, this refers to computing \code{nfolds} separate
statistics, and then using their mean and estimated standard
error to describe the CV curve. If \code{grouped=FALSE}, an error
matrix is built up at the observation level from the predictions
from the \code{nfold} fits, and then summarized (does not apply to
\code{type.measure="auc"}). For the \code{"cox"} family, \code{grouped=TRUE}
obtains the CV partial likelihood for the Kth fold by
\emph{subtraction}; by subtracting the log partial likelihood
evaluated on the full dataset from that evaluated on the on the
(K-1)/K dataset. This makes more efficient use of risk sets. With
\code{grouped=FALSE} the log partial likelihood is computed only on
the Kth fold}

\item{keep}{If \code{keep=TRUE}, a \emph{prevalidated} array is returned
containing fitted values for each observation and each value of
\code{lambda}. This means these fits are computed with this
observation and the rest of its fold omitted. The \code{foldid} vector
is also returned.  Default is keep=FALSE.}

\item{trace.it}{If \code{trace.it=1}, then progress bars are displayed;
useful for big models that take a long time to fit.}

\item{\dots}{Other arguments that can be passed to \code{multiview}}
}
\value{
an object of class \code{"cv.multiview"} is returned, which is a
list with the ingredients of the cross-validation
fit. \item{lambda}{the values of \code{lambda} used in the fits.}
\item{cvm}{The mean cross-validated error - a vector of length
\code{length(lambda)}.} \item{cvsd}{estimate of standard error of
\code{cvm}.} \item{cvup}{upper curve = \code{cvm+cvsd}.} \item{cvlo}{lower
curve = \code{cvm-cvsd}.} \item{nzero}{number of non-zero coefficients
at each \code{lambda}.} \item{name}{a text string indicating type of
measure (for plotting purposes).} \item{multiview.fit}{a fitted
multiview object for the full data.} \item{lambda.min}{value of
\code{lambda} that gives minimum \code{cvm}.} \item{lambda.1se}{largest
value of \code{lambda} such that error is within 1 standard error of
the minimum.} \item{fit.preval}{if \code{keep=TRUE}, this is the array
of prevalidated fits. Some entries can be \code{NA}, if that and
subsequent values of \code{lambda} are not reached for that fold}
\item{foldid}{if \code{keep=TRUE}, the fold assignments used}
\item{index}{a one column matrix with the indices of \code{lambda.min}
and \code{lambda.1se} in the sequence of coefficients, fits etc.}
}
\description{
Does k-fold cross-validation (CV) for multiview and produces a CV curve.
}
\details{
The current code can be slow for "large" data sets, e.g. when the
number of features is larger than 1000.  It can be helpful to see
the progress of multiview as it runs; to do this, set trace.it = 1
in the call to multiview or cv.multiview.  With this, multiview
prints out its progress along the way.  One can also pre-filter the
features to a smaller set, using the exclude option, with a filter
function.

If there are missing values in the feature matrices:
we recommend that you center the columns of each feature matrix, and then fill in the missing values with 0.

For example, \cr
\code{x <- scale(x,TRUE,FALSE)} \cr
\code{x[is.na(x)] <- 0} \cr
\code{z <- scale(z,TRUE,FALSE)} \cr
\code{z[is.na(z)] <- 0} \cr

Then run multiview in the usual way. It will exploit the assumed shared latent factors
to make efficient use of the available data.

The function runs \code{multiview} \code{nfolds+1} times; the first to get the
\code{lambda} sequence, and then the remainder to compute the fit with each
of the folds omitted. The error is accumulated, and the average error and
standard deviation over the folds is computed.  Note that \code{cv.multiview}
does NOT search for values for \code{rho}. A specific value should be
supplied, else \code{rho=0} is assumed by default. If users would like to
cross-validate \code{rho} as well, they should call \code{cv.multiview} with
a pre-computed vector \code{foldid}, and then use this same fold vector in
separate calls to \code{cv.multiview} with different values of \code{rho}.
}
\examples{
# Gaussian
# Generate data based on a factor model
set.seed(1)
x = matrix(rnorm(100*20), 100, 20)
z = matrix(rnorm(100*20), 100, 20)
U = matrix(rnorm(100*5), 100, 5)
for (m in seq(5)){
    u = rnorm(100)
    x[, m] = x[, m] + u
    z[, m] = z[, m] + u
    U[, m] = U[, m] + u}
x = scale(x, center = TRUE, scale = FALSE)
z = scale(z, center = TRUE, scale = FALSE)
beta_U = c(rep(0.1, 5))
y = U \%*\% beta_U + 0.1 * rnorm(100)
fit1 = cv.multiview(list(x=x,z=z), y, rho = 0.3)

# plot the cross-validation curve
plot(fit1)

# extract coefficients
coef(fit1, s="lambda.min")

# extract ordered coefficients
coef_ordered(fit1, s="lambda.min")

# make predictions
predict(fit1, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min")

# Binomial
\donttest{
by = 1 * (y > median(y)) 
fit2 = cv.multiview(list(x=x,z=z), by, family = binomial(), rho = 0.9)
predict(fit2, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min", type = "response")
plot(fit2)
coef(fit2, s="lambda.min")
coef_ordered(fit2, s="lambda.min")

# Poisson
py = matrix(rpois(100, exp(y))) 
fit3 = cv.multiview(list(x=x,z=z), py, family = poisson(), rho = 0.6)
predict(fit3, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min", type = "response") 
plot(fit3)
coef(fit3, s="lambda.min")
coef_ordered(fit3, s="lambda.min")
}
}
