#' @method family multiview
#' @export
family.multiview  <- function(object, ...) {
    ## families=c(elnet = "gaussian", lognet = "binomial", fishnet = "poisson",
    ##            multnet = "multinomial", coxnet = "cox", mrelnet = "mgaussian")
    ## cl <- class(object)[1]
    ## families[cl]
    fam <- object$family
    if (is.character(fam)) fam else fam$family
}

#' @method family cv.multiview
#' @export
family.cv.multiview <- function(object, ...) family(object$multiview.fit)
