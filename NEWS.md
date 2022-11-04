# CVXR 1.0-11

* Being more careful about coercing to `dgCMatrix`
  via `(as(as(<matrix>, "CsparseMatrix"), "generalMatrix"))`
* Modify all class inheritance checks to use `inherits()`

# CVXR 1.0-10

* Now requiring the updated scs 3.0 as an import

# multiview 0.7
  * Added Cox model
  * Like `glmnet` we now use `standardize = TRUE` by default.
  
# multiview 0.5 to 0.6
  * Internal development releases (not on CRAN)
  
# multiview 0.4

* First CRAN release (2022-09-02).
