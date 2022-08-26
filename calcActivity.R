#' A function to calculate activity level.
#' Activity id defined as the euclidean distance of a sample from 0.
#'
#' @param assay A matrix of normalized feature values in a long-format.
#' Columns of assay must include 'experiment', 'plate', 'well', 'MOA',
#' 'compound' , 'concentration' (uM) and 'value' for type=CP0.
#' For type=LC, columns of assay must include 'experiment', 'plate', 'well'  and 'value'
#' @param type the type of the plate for the assay. Possible values are 'CP70' for CP70 plates
#' and "LC" for low controls plates
#'
#' @return a data frame
#'
calcActivity <- function(assay, type){
  if(type=="CP70"){
    activity <- aggregate(assay$value,
                             list(experiment=assay$experiment,
                                  plate=assay$plate,
                                  well=assay$well,
                                  MOA=assay$MOA,
                                  compound=assay$compound,
                                  concentration=assay$concentration),
                             FUN=function(y) sqrt(sum(y^2)))
  }else if(type=="LC"){
    activity <- aggregate(assay$value,
                            list(experiment=assay$experiment,
                                 plate=assay$plate,
                                 well=assay$well),
                            FUN=function(y) sqrt(sum(y^2)))
    activity$MOA=NA
    activity$compound=NA
    activity$concentration = 0
    activity <- activity[, c("experiment", "plate", "well", "MOA",
                             "compound", "concentration", "x")]
  }
  activity$type <- type
  return(activity)
}
