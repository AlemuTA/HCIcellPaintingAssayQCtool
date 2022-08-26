#' Fisher Z transformation of a correlation r
#'
#' @param r a numeric value (correlation coefficient)
#'
#' @return a numeric value

fishZtransform <- function(r){
  0.5*log10((1+r)/(1-r))
  }
