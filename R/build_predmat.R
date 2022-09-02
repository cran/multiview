build_predmat <- function(outlist, lambda, x_list, offset, foldid, alignment, family, ...){
  ## Watch out: glmnet dispatches on methods for this function...
  m  <- length(x_list)
  N  <- nrow(x_list[[1L]])
  if (!is.null(offset)) {
    is.offset <- TRUE
    offset <- drop(offset)
  } else {
    is.offset = FALSE
  }

  predmat <- matrix(NA, N, length(lambda))
  nfolds <- max(foldid)
  nlams <- double(nfolds)
  nlambda <- length(lambda)

  for (i in seq(nfolds)) {
    which <- (foldid == i)
    fitobj = outlist[[i]]
    x_sub_list <- lapply(x_list, function(x) x[which, , drop = FALSE])
    if (is.null(offset)) {
      offset_sub  <- NULL
    } else {
      offset_sub <- rep(offset_sub[which, drop = FALSE])
    }
    preds <- switch(alignment,
                    fraction = predict(fitobj, x_sub_list, newoffset = offset_sub,...),
                    lambda = predict(fitobj, x_sub_list, s = lambda, newoffset = offset_sub,...)
                    )
    nlami <- min(ncol(preds), nlambda)
    predmat[which, seq(nlami)]  <- preds[, seq(nlami)]
    if (nlami < nlambda) {
      predmat[which, seq(from = nlami, to = nlambda)] <- preds[, nlami]
    }
  }
  rn <- rownames(x_list[[1L]])
  sn <- paste("s", seq(0, length = nlambda), sep = "")
  dimnames(predmat) <- list(rn, sn)
  attr(predmat, "family") <- family
  predmat
}

