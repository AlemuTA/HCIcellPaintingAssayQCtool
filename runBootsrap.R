#' A function to run parametric-bootsrap
#'
#' @param train.assay a data frame for the training CP70 plates (output from the function 
#' \code{prepareInputRawData()})

#' @return a data frame for the simulated plates
runBootsrap <- function(train.assay){
  train.assay.CP70 <- train.assay$CP70.Z
  train.assay.summarized <- aggregate(train.assay.CP70$value,
                                     by=list(experiment=train.assay.CP70$experiment,
                                             plate=train.assay.CP70$plate,
                                             compound_cnctr=train.assay.CP70$compound_cnctr,
                                             feature=train.assay.CP70$feature),
                                     FUN=mean, na.rm=TRUE)


  train.assay.summarized <-
    train.assay.summarized[train.assay.summarized$compound_cnctr %in% selecteCpd$compound_cnctr,]



  train.assay.summarized.melt <- reshape2::dcast(train.assay.summarized, formula =
                              plate+experiment+compound_cnctr ~ feature, value.var="x")

  bootsrap_dat <- do.call(rbind,
        lapply(unique(train.assay.summarized.melt$compound_cnctr),
               function(moa_cpd){
                 message(paste0("... ",moa_cpd))

                 Y <- train.assay.summarized.melt[train.assay.summarized.melt$compound_cnctr==moa_cpd,]
                 #saveRDS(Y, paste0("temp/data_", gsub("/","-", moa_cpd), ".rds"))
                 YB <- individualPlateBootsrap(Y)
                 YB
               }))
  bootsrap_dat
}
