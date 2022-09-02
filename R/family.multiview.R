#' @method family multiview
#' @export
family.multiview  <- function(object, ...) {
    families=c(elnet = "gaussian", lognet = "binomial", fishnet = "poisson",
               multnet = "multinomial", coxnet = "cox", mrelnet = "mgaussian")
    cl <- class(object)[1]
    families[cl]
}

#' @method family cv.multiview
#' @export
family.cv.multiview <- function(object, ...) family(object$glmnet.fit)
