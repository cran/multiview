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

#' Collapse a list of named lists into one list with the same name
#' @param in_list a list of named lists all with same names (not checked for efficiency)
#' @return a single list with named components all concatenated
collapse_named_lists  <- function(in_list) {
  lnames  <- names(in_list[[1L]])
  result  <- lapply(lapply(lnames, function(x) lapply(in_list, function(y) y[[x]])), unlist)
  names(result)  <- lnames
  result
}


#' Translate indices in `1:nvars` to column indices in list of x
#' matrices. No sanity checks
#' @inheritParams multiview
#' @param index vector of indices between 1 and nvars = sum of
#'   `ncol(x)` for x in x_list
#' @return a conformed list of column indices for each matrix,
#'   including possibly column indices of length 0
# @examples
# x_list  <- list(x = matrix(1:12, nrow = 3), y = matrix(1:9, nrow = 3),
#                 z = matrix(1:24, nrow = 3), w = matrix(1:30, nrow = 3))
# i  <- sort(sample(25L, size = 5L))
# multiview:::to_xlist_index(x_list, i)
#
to_xlist_index <- function(x_list, index) {
  m  <- length(x_list)
  p_x <- sapply(x_list, ncol)
  starts  <- cumsum(c(1L, p_x[-m]))
  p_index  <- mapply(seq.int, starts, length.out = p_x, SIMPLIFY = FALSE)
  x_list_indices <- lapply(p_index, intersect, index)
  grouped_indices <- mapply(function(x, y) x - y, x_list_indices, starts - 1L, SIMPLIFY = FALSE)
  names(grouped_indices)  <- names(x_list)
  grouped_indices
}

#' Translate from column indices in list of x matrices to indices in
#' `1:nvars`. No sanity checks for efficiency
#' @inheritParams multiview
#' @param index_list a list of column indices for each matrix,
#'   including possibly column indices of length 0
#' @return a vector of indices between 1 and nvars = sum of `ncol(x)`
#'   for x in x_list
# @examples
# x_list  <- list(x = matrix(1:12, nrow = 3), y = matrix(1:9, nrow = 3),
#                 z = matrix(1:24, nrow = 3), w = matrix(1:30, nrow = 3))
# i  <- sort(sample(25, size = 10))
# ix <- multiview:::to_xlist_index(x_list, i)
# j  <- multiview:::to_nvar_index(x_list, ix)
# all(i == unlist(j))  # true
# all(mapply(function(x, y) all(x == y), ix, multiview:::to_xlist_index(x_list, unlist(j)))) # true
#
to_nvar_index <- function(x_list, index_list) {
  m  <- length(x_list)
  stopifnot(m == length(index_list))
  p_x <- sapply(x_list, ncol)
  ends  <- cumsum(c(0L, p_x[-m]))
  mapply(function(x, y) x + y, index_list, ends, SIMPLIFY = FALSE)
}


#' Select x_list columns specified by (conformable) list of indices
#' @inheritParams multiview
#' @param indices a vector of indices in `1:nvars`
#' @return a list of x matrices
# @examples
# x_list  <- list(x = matrix(1:12, nrow = 3), y = matrix(1:9, nrow = 3),
#                 z = matrix(1:24, nrow = 3), w = matrix(1:30, nrow = 3))
# multiview:::select_matrix_list_columns(x_list, list(1:2, 1:2, 3:4, 6:10))
select_matrix_list_columns  <- function(x_list, indices) {
  stopifnot(length(x_list) == length(indices))
  mapply(function(x, y) x[, y], x_list, indices, SIMPLIFY = FALSE)
}

#' Return a new list of x matrices of same shapes as those in x_list
#' @inheritParams multiview.fit
reshape_x_to_xlist <- function(x, x_list) {
  m  <- length(x_list)
  p_x <- sapply(x_list, ncol)
  ends  <- cumsum(c(0L, p_x[-m]))
  indices  <- mapply(function(p, pad) seq_len(p) + pad, p_x, ends, SIMPLIFY = FALSE)
  lapply(indices, function(j) x[, j])
}

