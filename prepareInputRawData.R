#' A function to prepare input raw data (test or training data)
#'
#' @param list.datapath A list of paths where to extract the raw data. The source data must be in CSV format and the
#' first few columns (annotation columns) must include "PLATE_ID", "EXPERIMENT_NAME", "WELL_ID", "SALT_FORM",
#' and "CONCENTRATION", and the last columns must be the features
#' @param MOA.features a data frame containing the list of selected features
#' @param cpdAnnot a data frame of annotation for the compounds (SALTFORM) and mode-of-action (MOA)
#'
#' @return a list of data frames: one for CP70 another one for the low-control (LC) samples
prepareInputRawData <- function(list.datapath, MOA.features, cpdAnnot){
  #
  assay.all.plates <- lapply(list.datapath$datapath, function(plt){
    print(plt)
    # read experiment
    temp <- read.csv(plt)
    required.columns <- c("PLATE_ID", "EXPERIMENT_NAME", "WELL_ID",
                          "SALT_FORM", "CONCENTRATION")
    # subset samples and the MOA features
    ## CP70
    CP70 <- temp[temp$WELLTYPE_CODE=="SAMPLE",
                 c(required.columns, MOA.features$feature)]

    ## LC
    LC <- temp[temp$WELLTYPE_CODE=="LC",
               c(required.columns, MOA.features$feature)]

    # # remove invalidated wells
    # CP70 <- CP70[CP70$IS_VALID %in% c(0,1), ]
    # LC   <- LC[LC$IS_VALID %in% c(0,1), ]


    # standardize the samples
    ## CP70
    CP70.Z <- sapply(colnames(CP70)[-c(1:5)], function(x){
      y <- (CP70[, x] - mean(LC[, x], na.rm = TRUE))/sd(LC[, x], na.rm = TRUE)
      y
    })
    ## LC
    LC.Z <- sapply(colnames(LC)[-c(1:5)], function(x){
      y <- (LC[, x] - mean(LC[, x], na.rm = TRUE))/sd(LC[, x], na.rm = TRUE)
      y
    })
    # combine annotations and standarzed samples
    CP70.Z <- cbind(CP70[, 1:5], CP70.Z)
    LC.Z   <- cbind(LC[, 1:5], LC.Z)

    # summarize across replicates (plates and wells)
    ## reformat to a long form
    CP70.Z.melt <- reshape2::melt(CP70.Z, id.vars=required.columns)
    LC.Z.melt <- reshape2::melt(LC.Z, id.vars=required.columns)
    # add MOA info
    colnames(CP70.Z.melt) <- colnames(LC.Z.melt) <-
      c("plate", "experiment", "well", "compound", "concentration", "feature", "value")
    CP70.Z.melt$MOA <- sapply(CP70.Z.melt$compound, function(cpd){
      y <- cpdAnnot$MoA.Short[cpdAnnot$JNJ.Salt==cpd]
      if(length(y)==0){
        NA
      }else if(length(unique(y))>1){
        "conflicting!"
      }else{
        unique(y)
      }
    })
    CP70.Z.melt <- CP70.Z.melt[, c("plate", "experiment", "well",
                                   "MOA", "compound", "concentration",
                                   "feature", "value")]

    LC.Z.melt$MOA <- LC.Z.melt$compound<- LC.Z.melt$concentration <-  NA
    LC.Z.melt <- LC.Z.melt[, c("plate", "experiment", "well",
                                   "MOA", "compound", "concentration",
                                   "feature", "value")]

    CP70.Z.melt$concentration <- signif(CP70.Z.melt$concentration*1e6, digits = 3)

    # return output
    list(CP70=CP70.Z.melt, LC=LC.Z.melt)
  })

  CP70.Z <- do.call(rbind, lapply(assay.all.plates, function(x) x$CP70))
  LC.Z   <- do.call(rbind, lapply(assay.all.plates, function(x) x$LC))

  CP70.Z <- CP70.Z[!is.na(CP70.Z$MOA), ]
  CP70.Z$MOA <- sapply(CP70.Z$MOA, as.character)
  CP70.Z$compound_cnctr <-
    paste0("(",CP70.Z$MOA, ") ",  CP70.Z$compound, "_@", CP70.Z$concentration, "uM")

  # CP70.Z$experiment <- sapply(CP70.Z$experiment,
  #                             function(x) unlist(strsplit(x, "1536"))[2])
  # LC.Z$experiment <- sapply(LC.Z$experiment,
  #                           function(x) unlist(strsplit(x, "1536"))[2])
  list(CP70.Z=CP70.Z, LC.Z=LC.Z)
}
