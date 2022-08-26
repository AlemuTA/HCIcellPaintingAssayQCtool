#' A function to calculate the reference biosignature
#'
#' @param bootsrapData a data frame (output from \code{runBootsrap()})
#' @param summary.stat a name of function to calculate the reference signature. By default, summary.stat=mean, i.e.
#' the reference signature for a particular treatment is the mean signature across simulated plates
#'
#' @return a data frame
calcRefBiosig <- function(bootsrapData, summary.stat=mean){
  ref.biosignature <- aggregate(bootsrapData$value,
                                by = list("compound_cnctr"=bootsrapData$compound,
                                          "feature"=bootsrapData$feature),
                                          FUN=match.fun(summary.stat))
  ref.biosignature
}
