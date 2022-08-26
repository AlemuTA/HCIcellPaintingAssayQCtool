#' A wrapper function to process the training experiments
#'
#' @param NewTrainingPlates A data.frame that contains one row for each selected file, and following columns:
#' name: The filename. This is not the path to read to get at the actual data that was uploaded.
#' size: The size of the uploaded data, in bytes.
#' type: The MIME type r(for example, text/plain), or empty string.
#' datapath: The path to a temp file that contains the data that was uploaded.
#' @param cpdConcuM a CSV filecontaining the selected compounds and concentrations. Columns must include:
#' 'MOA', 'compound', 'Conc_uM1', and 'Conc_uM2'. The concentrations must be in micro-molar (uM).
#' @param MOA.features a CSV file with two columns: 'number' and 'feature'. Thefirst column is optional and it has no use
#' @param cpdAnnot a CSV file with columns at least 'JNJ Salt',	and	'MoA.Short'.
#' @param trainSaveFileName the name of the file to save the training results
#'
#' @return a list of results from processing the training plates.The function first saves the result
#' with a file name 'trainSaveFileName.rds' inside a sub folder named (temporarily) 'trainingRes'.
#'
processTrainingPlates <- function(NewTrainingPlates, cpdConcuM, MOA.features,
                                  cpdAnnot, trainSaveFileName){
  # read and prepare slected MOA features
  message("Reading inpputs  ...")
  MOA.features <- read.csv(MOA.features$datapath)
  MOA.features$feature <- gsub(".Zscore", "", MOA.features$feature)

  # read the CP70 plates annotation
  cpdAnnot <- read.csv(cpdAnnot$datapath)

  # calculate compound activity
  message("Preparing raw plate data ...")
  train.assay.all      <- prepareInputRawData(list.datapath = NewTrainingPlates, 
                                              MOA.features, cpdAnnot)
  message("Calculating compound acitivity ...")
  train.activity.CP70  <- calcActivity(assay = train.assay.all$CP70.Z, type = "CP70")
  train.activity.LC    <- calcActivity(assay = train.assay.all$LC.Z, type = "LC")
  train.activity       <- rbind(train.activity.CP70, train.activity.LC)

  # read and prepare the selcted compounds and concentration
  selecteCpd      <- read.csv(cpdConcuM$datapath)
  selecteCpd      <- processSelectedCpd(selecteCpd)

  # run bootsrap
  message("Running bootsrap ... (this takes a while)")
  bootsrapData    <- runBootsrap(train.assay = train.assay.all, selecteCpd = selecteCpd)

  # calculate reference biosignature
  message("Calculating reference signature ...")
  ref.biosignature<- calcRefBiosig(bootsrapData)

  # calculate QC scores per compound
  message("Calculating QC metrics ...")
  boot.QC.scores  <- calcQC(bootsrapData, selecteCpd, ref.biosignature)

  # summarize QC scores to a plate level and determine 95% QC limits
  plateQCprop     <- calcPlateQCprop(boot.QC.scores, boot.QC.scores)
  var.mean <- function(y){var(y)/length(y)}
  bootPlate.insideQClimit.prop <-
    list(zonePropMean=tapply(plateQCprop$value, plateQCprop$variable, mean),
         zonePropVar=tapply(plateQCprop$value, plateQCprop$variable, var.mean))

  # organize outputs, save and return
  train.results <- list(train.assay.all=train.assay.all,
                        bootsrapData = bootsrapData,
                        ref.biosignature = ref.biosignature,
                        boot.QC.scores = boot.QC.scores,
                        bootPlate.insideQClimit.prop = bootPlate.insideQClimit.prop,
                        train.activity = train.activity,
                        selecteCpd=selecteCpd,
                        MOA.features=MOA.features,
                        cpdAnnot=cpdAnnot)
  saveRDS(train.results, paste0("trainingRes/trainingRes_", trainSaveFileName, ".rds"))
  train.results
}
