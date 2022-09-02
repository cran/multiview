#' Make predictions from a "cv.multiview" object.
#'
#' This function makes predictions from a cross-validated multiview model, using
#' the stored `"multiview"` object, and the optimal value chosen for
#' `lambda`.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object Fitted `"cv.multiview"` or object.
#' @param newx List of new view matrices at which predictions are to be made.
#' @param s Value(s) of the penalty parameter `lambda` at which
#' predictions are required. Default is the value `s="lambda.1se"` stored
#' on the CV `object`. Alternatively `s="lambda.min"` can be used. If
#' `s` is numeric, it is taken as the value(s) of `lambda` to be
#' used. (For historical reasons we use the symbol 's' rather than 'lambda' to
#' reference this parameter)
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on the \dots{} argument which is passed
#' on to the `predict` method for `multiview` objects.
#'
#' @examples
#' # Gaussian
#' # Generate data based on a factor model
#' set.seed(1)
#' x = matrix(rnorm(100*10), 100, 10)
#' z = matrix(rnorm(100*10), 100, 10)
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
#' predict(fit1, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min")
#' 
#' # Binomial
#' \donttest{
#' by = 1 * (y > median(y)) 
#' fit2 = cv.multiview(list(x=x,z=z), by, family = binomial(), rho = 0.9)
#' predict(fit2, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min", type = "response")
#' 
#' # Poisson
#' py = matrix(rpois(100, exp(y))) 
#' fit3 = cv.multiview(list(x=x,z=z), py, family = poisson(), rho = 0.6)
#' predict(fit3, newx = list(x[1:5, ],z[1:5,]), s = "lambda.min", type = "response")
#' }
#' @importFrom stats predict
#' @method predict cv.multiview
#' @export
predict.cv.multiview = function(object, newx, s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s))
    lambda = s
  else
    if(is.character(s)) {
      s = match.arg(s)
      lambda = object[[s]]
      names(lambda) = s
    }
  else stop("Invalid form for s")
  predict(object$multiview.fit, newx, s =lambda, ...)
}
