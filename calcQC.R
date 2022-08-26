
#' A function that calculates the MOA QC
#'
#' @param test.assay a dataframe (assay). columns must include 'plate', 'experiment', 'MOA',
#' 'compound', 'concentration' (uM), 'feature' and 'value'. This is an output from \code{runBootsrap()} function.
#' @param selecteCpd a dataframe of the selected compounds and concentrations.
#' Columns must be 'MOA', 'compound' and 'concentration' (uM). This is an output from  \code{processSelectedCpd()} function.
#' @param ref.sig a data frame for reference biosignature determined from the bootstrap
#' result using the training plates. Columns must be 'MOA', 'compound', 'concentration' (uM),
#' 'feature' and 'x'. (x is for the value column). This an output from \code{calcRefBiosig()} function.
#'
#' @return  a data frame
#'
calcQC <- function(test.assay, selecteCpd, ref.sig){
  test.assay.summarized <- aggregate(test.assay$value,
                       by=list(plate=test.assay$plate,
                               experiment=test.assay$experiment,
                               compound_cnctr=test.assay$compound_cnctr,
                               feature=test.assay$feature),
                       FUN=mean, na.rm=TRUE)


  test.assay.summarized <-
    test.assay.summarized[test.assay.summarized$compound_cnctr %in%
                            selecteCpd$compound_cnctr,]


  QC_scores <- do.call(rbind, lapply(unique(test.assay.summarized$compound_cnctr),
                                     function(moa_cpd){
        Yplat <- test.assay.summarized[test.assay.summarized$compound_cnctr == moa_cpd, ]
        Yref  <- ref.sig[ref.sig$compound_cnctr == moa_cpd,  ]

        scores <- as.data.frame(do.call(rbind,
              lapply(unique(Yplat$plate),  function(plt){
                x <- Yplat[Yplat$plate == plt,]
                xy <- merge(x, Yref,
                            by=c("compound_cnctr", "feature"))
                scr <- calcSimilarityMetrics(x=xy$x.x, y=xy$x.y)
                scr
              })))
        #scores[is.na(scores)]  <- 0
        scores$compound_cnctr <- moa_cpd
        scores$plate <- unique(Yplat$plate)
        scores$experiment <- sapply(scores$plate, function(plt){
          unique(Yplat$experiment[Yplat$plate==plt])
        })
        scores
  }))

  colnames(QC_scores)[1:2] <- c("Spearman rank corr.","Euclidian dist.")
  QC_scores$`Spearman rank corr.` <- fishZtransform(QC_scores$`Spearman rank corr.`)
  QC_scores$`Euclidian dist.` <- log(QC_scores$`Euclidian dist.`)
  QC_scores
}
