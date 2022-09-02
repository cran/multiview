#' Evaluate the contribution of data views in making prediction
#'
#' Evaluate the contribution of each data view in making prediction. The function has two options.
#' If `force` is set to `NULL`, the data view contribution is benchmarked by the null model.
#' If `force` is set to a list of data views, the contribution is benchmarked by the model fit on
#' this list of data views, and the function evaluates the marginal contribution of each additional data
#' view on top of this benchmarking list of views.
#' The function returns a table showing the percentage improvement in reducing error as compared to the bechmarking model
#' made by each data view.
#'
#' @inheritParams cv.multiview
#' @param eval_data If `train`, we evaluate the contribution of data views based on training data
#' using cross validation error; if `test`, we evaluate the contribution of data views based on test data.
#' Default is `train`. If set to `test`, users need to provide the test data, i.e.
#' `x_list_test` and `y_test`.
#' @param s Value(s) of the penalty parameter `lambda` at which
#' predictions are required. Default is the value `s="lambda.1se"` stored
#' on the CV `object`. Alternatively `s="lambda.min"` can be used. If
#' `s` is numeric, it is taken as the value(s) of `lambda` to be
#' used. (For historical reasons we use the symbol 's' rather than 'lambda' to
#' reference this parameter)
#' @param x_list_test A list of `x` matrices in the test data for evaluation.
#' @param test_y The quantitative response in the test data with length equal to the
#' number of rows in each `x` matrix of the test data.
#' @param force If `NULL`, the data view contribution is benchmarked by the null model.
#' If users want to benchmark by the model fit on a specified list of data views, `force` needs to
#' be set to this list of benchmarking data views, i.e. a list of `x` matrices. The function then
#' evaluates the marginal contribution of each additional data, i.e. the data views in `x_list` but not in
#' `force`, on top of the benchmarking views.
#' @return a data frame consisting of the view, error metric, and percentage improvement.
#' @examples
#' set.seed(3)
#' # Simulate data based on the factor model
#' x = matrix(rnorm(200*20), 200, 20)
#' z = matrix(rnorm(200*20), 200, 20)
#' w = matrix(rnorm(200*20), 200, 20)
#' U = matrix(rep(0, 200*10), 200, 10) # latent factors
#' for (m in seq(10)){
#'     u = rnorm(200)
#'     x[, m] = x[, m] + u
#'     z[, m] = z[, m] + u
#'     w[, m] = w[, m] + u
#'     U[, m] = U[, m] + u}
#' beta_U = c(rep(2, 5),rep(-2, 5))
#' y = U %*% beta_U + 3 * rnorm(100)
#'
#' # Split training and test sets
#' smp_size_train = floor(0.9 * nrow(x))
#' train_ind = sort(sample(seq_len(nrow(x)), size = smp_size_train))
#' test_ind = setdiff(seq_len(nrow(x)), train_ind)
#' train_X = scale(x[train_ind, ])
#' test_X = scale(x[test_ind, ])
#' train_Z <- scale(z[train_ind, ])
#' test_Z <- scale(z[test_ind, ])
#' train_W <- scale(w[train_ind, ])
#' test_W <- scale(w[test_ind, ])
#' train_y <- y[train_ind, ]
#' test_y <- y[test_ind, ]
#' foldid = sample(rep_len(1:10, dim(train_X)[1]))
#'
#' # Benchmarked by the null model:
#' rho = 0.3
#' view.contribution(x_list=list(x=train_X,z=train_Z), train_y, rho = rho,
#'                   eval_data = 'train', family = gaussian())
#' view.contribution(x_list=list(x=train_X,z=train_Z), train_y, rho = rho,
#'                   eval_data = 'test', family = gaussian(),
#'                   x_list_test=list(x=test_X,z=test_Z), test_y=test_y)
#'
#' # Force option -- benchmarked by the model train on a specified list of data views:
#' view.contribution(x_list=list(x=train_X,z=train_Z,w=train_W), train_y, rho = rho,
#'                   eval_data = 'train', family = gaussian(), force=list(x=train_X))
#'
#' @importFrom glmnet cv.glmnet
#' @export
#'
view.contribution = function(x_list, y, family = gaussian(),
                             rho, s = c("lambda.min", "lambda.1se"),
                             eval_data = c('train', 'test'),
                             weights = NULL, offset = NULL, mvlambda = NULL,
                             type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
                             x_list_test = NULL, test_y = NULL,
                             nfolds = 10, foldid = NULL,
                             force = NULL,
                             ...) {

  type.measure =  match.arg(type.measure)
  s  <- match.arg(s)
  m <- length(x_list)
  N <- nrow(x_list[[1L]])
  p <- do.call(sum, lapply(x_list, ncol))
  view_names  <- names(x_list)
  y <- drop(y)

  ## if (!is.null(lambda) && length(lambda) < 2)
  ##   stop("Need more than one value of lambda for cv.glmnet")

  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else {
    nfolds <- max(foldid)
  }
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }

  if (is.null(weights)) {
    weights  <- rep(1.0, N)
  }

  if (eval_data == 'test' && (is.null(x_list_test) | is.null(test_y))) {
    stop("please provide test data if you want to evaluate the view contribution based on test data")
  }

  if (is.null(force)){
    # Benchmarked by the null model
    if (eval_data == 'test'){
      err_null = mean((test_y - mean(y))^2)
    } else if (eval_data == 'train'){
      err_null_fold = sapply(X = seq(max(foldid)),
                             FUN = function(i) {
                               ind <- foldid == i
                               mean((mean(y[!ind]) - y[ind])^2)
                             })
      err_null = mean(err_null_fold)
    }

    # Cooperative learning, all views
    full_fit_no_pf = cv.multiview(x_list, y = y, rho = rho, family = family,
                                    weights = weights, type.measure = type.measure,
                                    foldid = foldid)
    if (eval_data == 'test'){
      yhat_test_no_pf = predict.cv.multiview(full_fit_no_pf, x_list_test, s=s)
      err_all_view = mean((yhat_test_no_pf - test_y)^2)
    } else if (eval_data == 'train'){
      err_all_view = min(full_fit_no_pf$cvm)
    }

    # Use only one data view
    err_each_view_list = c()
    err_view_name = c()
    for (i in seq(m)){
      data_view = x_list[[i]]
      fit_data_view = glmnet::cv.glmnet(data_view, y, standardize = F, foldid = foldid)

      if (eval_data == 'test'){
        data_view_test = x_list_test[[i]]
        data_view_pred = predict(fit_data_view, data_view_test, s=s)
        err_each_view = mean((data_view_pred - test_y)^2)
      } else if (eval_data == 'train'){
        err_each_view = min(fit_data_view$cvm)
      }
      err_each_view_list = c(err_each_view_list, err_each_view)
      err_view_name = c(err_view_name, view_names[i])
    }

    err_list = c(err_null, err_each_view_list, err_all_view)
    err_list_contribution = (err_null - err_list) / err_null
    df_res = data.frame(view = c("null", err_view_name, "cooperative (all)"),
                        error = err_list,
                        percentage_improvement = err_list_contribution*100)
  } else {
    x_list_base <- force
    base_model <- cv.multiview(x_list_base, y = y, rho = rho, family = family,
                              weights = weights, type.measure = type.measure,
                              foldid = foldid)

    if (eval_data == 'test'){
      x_list_base_test = x_list_test[names(x_list_test) %in% names(force)]
      yhat_base_model = predict(base_model, x_list_base_test, s=s)
      err_base_model = mean((yhat_base_model - test_y)^2)
    } else if (eval_data == 'train'){
      err_base_model = min(base_model$cvm)
    }

    x_list_additional <- x_list[!(names(x_list) %in% names(force))] #x_list[!(x_list %in% force)]
    additional_view_names <- names(x_list_additional)
    m_add <- length(x_list_additional)
    err_each_view_list <- c()
    err_view_name <- c()

    for (i in seq(m_add)){
      x_list_base_add = append(x_list_base, list(x_list_additional[[i]]))

      each_model = cv.multiview(x_list_base_add, y = y, rho = rho, family = family,
                                weights = weights, type.measure = type.measure,
                                foldid = foldid)

      if (eval_data == 'test'){
        x_list_additional_test = x_list_test[!(names(x_list_test) %in% names(force))]
        x_list_base_add_test = append(x_list_base_test, list(x_list_additional_test[[i]]))
        yhat_each_model = predict(each_model, x_list_base_add_test, s=s)
        err_each_model = mean((yhat_each_model - test_y)^2)
      } else if (eval_data == 'train'){
        err_each_model = min(each_model$cvm)
      }

      err_each_view_list = c(err_each_view_list, err_each_model)
      err_view_name = c(err_view_name, additional_view_names[i])
    }

  err_list = c(err_base_model, err_each_view_list)
  err_list_contribution = (err_base_model - err_list) / err_base_model
  df_res = data.frame(view = c("baseline", err_view_name),
                      error = err_list,
                      percentage_improvement = err_list_contribution*100)
  }

  return(df_res)
}
