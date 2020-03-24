#' Generate data correlated by groups
#'
#' Generates `m_indep * m_corr` covariates, with a grouped
#' correlation matrix with group size `m_corr`.
#'
#' @param n Number of observations in generated data
#' @param m_indep Number of independent covariates
#' @param m_corr Number of covariates correlated in group
#' @param corr_mat Correlation matix, has to be of dim m_corr. If missing
#'    it uses a random correlation matrix.
#' @param seed Integer to use as random seed. If missing, no seed is set.
#' @examples
#' get_grouped_corr_data(10, 3, 2, seed=7)
#' @export
get_grouped_corr_data <- function(n, m_indep, m_corr, corr_mat, seed) {
  if (!missing(seed)) {
    set.seed(seed)
  }
  if (missing(corr_mat)) {
  }

  X0 <- matrix(rnorm(n*m_indep), nrow=n)
  colnames(X0) <- paste0("x", 1:m_indep)
  for (idx in 1:m_indep) {
    x <- X0[, idx, drop=FALSE]
    X_corr <- get_corr_data(m_corr, x, corr_mat=corr_mat)
    colnames(X_corr) <- c(colnames(x), paste0(colnames(x), "_corr_", 1:(m_corr-1)))
    if (idx == 1) {
      X <- X_corr
    } else {
      X <- cbind(X, X_corr)
    }
  }
  return(X)
}

#' Generating a random positive-definite matrix
#'
#' Generating a random positive-definite matrix with user-specified
#' positive eigenvalues. If eigenvalues are not specified, they are
#' generated from a uniform distribution.
#'
#' Function taken from https://stat.ethz.ch/pipermail/r-help/2008-February/153708
#'
#' @param n Dimension of the matrix
#' @param ev eigenvalue vector
get_posdef_matrix <- function (n, ev = runif(n, 0, 10))
{
Z <- matrix(ncol=n, rnorm(n^2))
decomp <- qr(Z)
Q <- qr.Q(decomp)
R <- qr.R(decomp)
d <- diag(R)
ph <- d / abs(d)
O <- Q %*% diag(ph)
Z <- t(O) %*% diag(ev) %*% O
return(Z)
}

#' Generate correlated data
#'
#' @param p Number of correlated covariates
#' @param x Generate covariates correlated to `x`.
#' @param corr_mat Correlation matix, has to be of dim `p`. If missing
#'    it uses a random correlation matrix.
#' @param seed Integer to use as random seed. If missing, no seed is set.
#' @examples
#' get_corr_data(10, seed=7)
#' @export
get_corr_data <- function(p, x=rnorm(100), corr_mat, seed) {
  n <- length(x)
  if (!missing(seed)) {
    set.seed(seed)
  }

  x_extra <- scale(matrix( rnorm(n*(p-1)), ncol=(p-1)))
  x_all <- cbind(scale(x),x_extra)

  # find the current correlation matrix
  c1 <- var(x_all)

  # cholesky decomposition to get independence
  chol1 <- solve(chol(c1))

  newx <-  x_all %*% chol1

  # create new correlation structure (zeros can be replaced with other vals)
  if (missing(corr_mat)) {
    newc <- get_posdef_matrix(p)
  } else {
    newc <- corr_mat
  }

  chol2 <- chol(newc)

  finalx <- newx %*% chol2 * sd(x) + mean(x)
  return(finalx)
}
