#' A function to process (convert wide format to long format and creat treatment) the input for the selected compounds and concentration
#'
#' @param selecteCpd a data frame that contains the provided slecetd compounds and concentrations
#'
#' @return a data frame
#'
processSelectedCpd <- function(selecteCpd){
  selecteCpd      <- reshape2::melt(selecteCpd, id.vars=c("MOA", "compound"))
  colnames(selecteCpd)[colnames(selecteCpd)=="value"] <- "concentration"
  selecteCpd <- selecteCpd[,c("MOA", "compound", "concentration")]
  selecteCpd$compound_cnctr <-
    paste0("(",selecteCpd$MOA, ") ",
           selecteCpd$compound, "_@",  selecteCpd$concentration, "uM")
  selecteCpd
}
