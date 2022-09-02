## ----include=FALSE------------------------------------------------------------
# the code in this chunk enables us to truncate the print output for each
# chunk using the `out.lines` option
# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
        
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

## ---- message=FALSE, results='hide'-------------------------------------------
library(multiview)

## -----------------------------------------------------------------------------
quick_example <- readRDS(system.file("exdata", "example.RDS", package = "multiview"))
x <- quick_example$x
z <- quick_example$z
y <- quick_example$y

## -----------------------------------------------------------------------------
fit <- multiview(list(x,z), y, family=gaussian(), rho=0)

## -----------------------------------------------------------------------------
x <- scale(x,TRUE,FALSE)
x[is.na(x)] <- 0
z <- scale(z,TRUE,FALSE)
z[is.na(z)] <- 0

## -----------------------------------------------------------------------------
plot(fit)

## ----out.lines = 10-----------------------------------------------------------
print(fit)

## -----------------------------------------------------------------------------
set.seed(1)
nx <- matrix(rnorm(5 * 50), 5, 50)
nz <- matrix(rnorm(5 * 50), 5, 50)
predict(fit, newx = list(nx, nz), s = c(0.1, 0.05))

## ----out.lines = 10-----------------------------------------------------------
coef(fit, s = 0.1)

## ----out.lines = 10-----------------------------------------------------------
coef_ordered(fit, s=0.1)

## -----------------------------------------------------------------------------
cvfit <- cv.multiview(list(x,z), y, rho=0.1, family = gaussian(),
                      type.measure = "mse", nfolds = 20)

## -----------------------------------------------------------------------------
plot(cvfit)

## ----out.lines = 10-----------------------------------------------------------
cvfit$lambda.min
coef(cvfit, s = "lambda.min")

## -----------------------------------------------------------------------------
predict(cvfit, newx = list(x[1:5,],z[1:5,]), s = "lambda.min")

## -----------------------------------------------------------------------------
quick_example <- readRDS(system.file("exdata", "example_contribution.RDS",
                                     package = "multiview"))
rho <- 0.3
view.contribution(x_list=list(x=quick_example$x, z=quick_example$z), quick_example$y,
                  rho = rho, family = gaussian(), eval_data = 'train')

## -----------------------------------------------------------------------------
view.contribution(x_list=list(x=quick_example$x, z=quick_example$z, w=quick_example$w),
                  quick_example$y, rho = rho, eval_data = 'train', family = gaussian(),
                  force=list(x=quick_example$x))

## -----------------------------------------------------------------------------
view.contribution(x_list = list(x = quick_example$x, z = quick_example$z),
                  quick_example$y, rho = rho,
                  x_list_test = list(x = quick_example$test_X, z = quick_example$test_Z),
                  test_y = quick_example$test_y, family = gaussian(), eval_data = 'test')

## -----------------------------------------------------------------------------
uvar <- function(x, means = FALSE) {
  # if means = TRUE, the means and variances are returned, 
  # otherwise just the variances
  m <- colMeans(x)
  n <- nrow(x)
  x <- x - outer(rep(1,n),m)
  v <- colSums(x^2) / (n - 1)
  if (means) list(mean = m, var = v) else v
}

uvar_multiple <- function(x_list) lapply(x_list, function(x) uvar(x))

vfilter <- function(x_list, q = 0.3, ...) {
    v <- uvar_multiple(x_list)
    lapply(v, function(x) which(x < quantile(x, q)))
}

## -----------------------------------------------------------------------------
set.seed(1)
x <- matrix(rnorm(100 * 20), 100, 20)
z <- matrix(rnorm(100 * 20), 100, 20)
y <- rnorm(100)
fit.filter <- multiview(list(x,z), y, exclude = vfilter)

## -----------------------------------------------------------------------------
cvfit.filter <- cv.multiview(list(x,z), y, exclude = vfilter)

## ---- eval=FALSE--------------------------------------------------------------
#  fit <- multiview(list(x,z), y, trace.it = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  fit <- cv.multiview(list(x,z), y, trace.it = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  multiview.control(itrace = 1)

## -----------------------------------------------------------------------------
example_binomial <- readRDS(system.file("exdata", "example_binomial.RDS",
                                        package = "multiview"))
x <- example_binomial$x
z <- example_binomial$z
y <- example_binomial$by

## -----------------------------------------------------------------------------
multiview.control(mxitnr = 100)
fit <- multiview(list(x=x,z=z), y, family = binomial(), rho = 0)
multiview.control(factory = TRUE)

## -----------------------------------------------------------------------------
predict(fit, newx = list(x[5:10,],z[5:10,]), type = "class", s = c(0.05, 0.01))

## ---- warning=FALSE-----------------------------------------------------------
cvfit <- cv.multiview(list(x=x,z=z), y, family = binomial(),
                      type.measure = "deviance", rho = 0.5)

## -----------------------------------------------------------------------------
plot(cvfit)

## -----------------------------------------------------------------------------
example_poisson <- readRDS(system.file("exdata", "example_poisson.RDS",
                                       package = "multiview")) 
x <- example_poisson$x
z <- example_poisson$z
y <- example_poisson$py

## -----------------------------------------------------------------------------
fit <- multiview(list(x=x, z=z), y, family = poisson(), rho = 0)

## -----------------------------------------------------------------------------
plot(fit)

## ----out.lines = 7------------------------------------------------------------
coef(fit, s = 1)
predict(fit, newx = list(x[1:5,],z[1:5,]), type = "response", s = c(0.1,1))

## ---- warning = FALSE, error = TRUE-------------------------------------------
cvfit <- cv.multiview(list(x,z), y, family = poisson(), rho = 0.1)

## -----------------------------------------------------------------------------
multiview.control(mxitnr = 100)
cvfit <- cv.multiview(list(x,z), y, family = poisson(), rho = 0.1)

## -----------------------------------------------------------------------------
multiview.control(factory = TRUE)

