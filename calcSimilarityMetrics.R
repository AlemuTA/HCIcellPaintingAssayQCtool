#' A function that claculates the similarity metrics (spearman rank correlation
#' and euclidean distance)
#'
#' @param x a vnumeric vector of size n
#' @param y a vnumeric vector of size n
#'
#' @return a vector of numeric vector size 2
#'
#' @examples
calcSimilarityMetrics <- function(x, y){


  #spearman_cor <- wCorr::weightedCorr(x, y, weights = w, method = "Spearman")
  spearman_cor <- cor(x, y, method = "spearman")


  # Euclidean distance
  w <-  abs(y)/max(abs(y))
  euclidean_dist <- sqrt(sum(((x-y)^2)*w))

  c(spearman_cor=spearman_cor,euclidean_dist=euclidean_dist)
}
