build_predmat_coxnetlist <-
  function(outlist, lambda, x_list, offset, foldid, alignment, family, y,weights,grouped,type.measure="deviance",...){
    nfolds <- max(foldid)
    if ((length(weights)/nfolds < 10) && !grouped) grouped <- TRUE
    devtrue <- type.measure == "deviance"
    nlambda <- length(lambda)    
    cvraw <- if(devtrue) matrix(NA, nfolds, nlambda) else NULL
    
    m <- length(x_list)
    N <- nrow(x_list[[1L]])
    
    predmat <- matrix(NA, N, nlambda)
    rn <- rownames(x_list[[1L]])
    sn <- paste("s",seq(0,length=nlambda),sep="")
    dimnames(predmat) <- list(rn,sn)

    full_x <- do.call(cbind, x_list)
    for (i in seq(nfolds)) {
      which <- (foldid == i)
      fitobj <- outlist[[i]]
      coefmat <- switch(alignment,
                        fraction = predict(fitobj, type = "coefficients", ...),
                        lambda = predict(fitobj, type = "coefficients", s = lambda, ...)
                        )
      nlami <- min(ncol(coefmat), nlambda)

      if(devtrue) {
        if (grouped) {
          plfull <- glmnet::coxnet.deviance(x = full_x, y = y, offset = offset,
                                    weights = weights, beta = coefmat)
          plminusk <- glmnet::coxnet.deviance(x = full_x[!which, ],
                                      y = y[!which, ],
                                      offset = offset[!which],
                                      weights = weights[!which], 
                                      beta = coefmat)
          cvraw[i, seq(nlami)] <- (plfull - plminusk)[seq(nlami)]
        } else {
          plk <- glmnet::coxnet.deviance(x = full_x[which, ],
                                 y = y[which, ], 
                                 offset = offset[which], 
                                 weights = weights[which],
                                 beta = coefmat)
          cvraw[i, seq(nlami)] <- plk[seq(nlami)]
        }
      }
      predmat[which, seq(nlami)] <- as.matrix(full_x[which,] %*% coefmat)
      if (nlami < nlambda){
        if(devtrue) cvraw[i, seq(from = nlami, to = nlambda)] <- cvraw[i, nlami]
        predmat[which,seq(from=nlami,to=nlambda)] <- predmat[which,nlami]
      }
    }
    if(devtrue) attr(predmat,"cvraw") <- cvraw
    attr(predmat, "family") <- family
    predmat
  }
