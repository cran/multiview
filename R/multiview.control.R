#' Internal multiview parameters
#'
#' View and/or change the factory default parameters in multiview
#'
#' If called with no arguments, `multiview.control()` returns a list with the
#' current settings of these parameters. Any arguments included in the call
#' sets those parameters to the new values, and then silently returns. The
#' values set are persistent for the duration of the R session.
#'
#' @param fdev minimum fractional change in deviance for stopping path; factory
#' default = 1.0e-5
#' @param devmax maximum fraction of explained deviance for stopping path;
#' factory default = 0.999
#' @param eps minimum value of lambda.min.ratio (see multiview); factory default=
#' 1.0e-6
#' @param big large floating point number; factory default = 9.9e35. Inf in
#' definition of upper.limit is set to big
#' @param mnlam minimum number of path points (lambda values) allowed; factory
#' default = 5
#' @param pmin minimum probability for any class. factory default = 1.0e-9.
#' Note that this implies a pmax of 1-pmin.
#' @param exmx maximum allowed exponent. factory default = 250.0
#' @param prec convergence threshold for multi response bounds adjustment
#' solution. factory default = 1.0e-10
#' @param mxit maximum iterations for multiresponse bounds adjustment solution.
#' factory default = 100
#' @param itrace If 1 then progress bar is displayed when running `multiview`
#' and `cv.multiview`. factory default = 0
#' @param epsnr convergence threshold for `multiview.fit`. factory default =
#' 1.0e-6
#' @param mxitnr maximum iterations for the IRLS loop in `multiview.fit`. factory
#' default = 25
#' @param factory If `TRUE`, reset all the parameters to the factory
#' default; default is `FALSE`
#' @return A list with named elements as in the argument list
#' @seealso `multiview`
#' @keywords models regression
#' @importFrom glmnet glmnet.control
#' @examples
#'
#' multiview.control(fdev = 0)  #continue along path even though not much changes
#' multiview.control()  # view current settings
#' multiview.control(factory = TRUE)  # reset all the parameters to their default
#'
#' @export multiview.control
multiview.control <- glmnet::glmnet.control
