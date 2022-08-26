# The following pacakges are taken from the SIBER package
pointsToEllipsoidCustom <- function(X, Sigma, mu) {
  if (ncol(Sigma) != nrow(Sigma))
    stop("Sigma must be a square matrix")
  if (ncol(X) != ncol(Sigma))
    stop("number of columns in X must be of same dimension as Sigma")
  if (length(mu) != ncol(Sigma))
    stop("length of mu must be of same dimension as Sigma")
  eig <- eigen(Sigma)
  SigSqrt = eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Z <- t(apply(X, 1, ellipsoidTransformCustom, SigSqrt, mu))
  return(Z)
}

ellipsoidTransformCustom <- function (x, SigSqrt, mu) {
  return(solve(SigSqrt, x - mu))
}


ellipseInOutCustom <- function (Z, p = 0.95, r = NULL)
{
  if (is.null(r)) {
    r <- stats::qchisq(p, df = ncol(Z))
  }
  inside <- rowSums(Z^2) < r
  return(inside)
}
