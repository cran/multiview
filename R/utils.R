#' Build a block row matrix for multiview
#' @param x_list list of x matrices
#' @param p_x a list of ncol of elements in x_list
#' @param pair an integer vector of two indices
#' @param rho the rho value
#' @return a block row of matrix for multiview
make_row  <- function(x_list, p_x, pair, rho) {
  sqrt_rho <- sqrt(rho); i <- pair[1]; j <- pair[2];
  result  <- lapply(p_x, function(p) matrix(0, nrow = nrow(x_list[[1L]]), ncol = p))
  result[[i]] <- -sqrt_rho * x_list[[i]]
  result[[j]] <- sqrt_rho * x_list[[j]]
  do.call(cbind, result)
}

