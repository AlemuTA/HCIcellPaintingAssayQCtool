#' A function to calculate plate-level QC metrics (a summary across treatments).
#'
#' @param QC_scores a data frame. Columns must include: 'Spearman rank corr.', 'Euclidian dist.',
#'  'compound_cnctr', 'plate' and 'experiment'. This is an output from \code{calcQC()} function.
#' @param boot_QC_scores a data frame. Columns must include: 'Spearman rank corr.', 'Euclidian dist.', 'plate', and
#' 'compound_cnctr. This is an output from \code{calcQC()} function.
#'
#' @return a data frame
calcPlateQCprop <- function(QC_scores, boot_QC_scores, sig.levl=0.05){

  sig.level.adj <- 1-sig.levl/length(unique(QC_scores$compound_cnctr)) # Benferroni

  plate.insideQClimit <-
    do.call(rbind,
           lapply(unique(QC_scores$compound_cnctr), function(moa_cpd){
             Yb <- boot_QC_scores[boot_QC_scores$compound_cnctr==moa_cpd,
                                  c("Spearman rank corr.", "Euclidian dist.",
                                    "plate", "compound_cnctr")]
             Yo <- QC_scores[QC_scores$compound_cnctr==moa_cpd,
                             c("Spearman rank corr.", "Euclidian dist.",
                               "plate", "experiment","compound_cnctr")]
             Z <- pointsToEllipsoidCustom(Yo[, 1:2], var(Yb[, 1:2]),
                                          colMeans(Yb[, 1:2]))
             data.frame(experiment=Yo$experiment,
                        plate=Yo$plate,
                        compound_cnctr = moa_cpd,
                        red_zone=!ellipseInOutCustom(Z, p =sig.level.adj) &
                          (Yo[, 1]< median(Yb[, "Spearman rank corr."])) &
                          (Yo[, 2]> median(Yb[, "Euclidian dist."])),

                        orange_zone=!ellipseInOutCustom(Z, p =sig.level.adj) &
                          (((Yo[, 1]> median(Yb[, "Spearman rank corr."])) &
                              (Yo[, 2]> median(Yb[, "Euclidian dist."]))) |
                             ((Yo[, 1]< median(Yb[, "Spearman rank corr."])) &
                                (Yo[, 2]< median(Yb[, "Euclidian dist."])))),

                        green_zone=ellipseInOutCustom(Z, p =sig.level.adj) |
                          ((Yo[, 1]>= median(Yb[, "Spearman rank corr."])) &
                             (Yo[, 2]<= median(Yb[, "Euclidian dist."]))))
           }))
  plate.insideQClimit.prop <-
    aggregate(plate.insideQClimit[, c("red_zone", "orange_zone", "green_zone")],
              by=list(plate=plate.insideQClimit$plate,
                      experiment = plate.insideQClimit$experiment),
              FUN=mean)

  plate.insideQClimit.prop <-
    plate.insideQClimit.prop[order(plate.insideQClimit.prop$green_zone,                                               plate.insideQClimit.prop$orange_zone,                                                    plate.insideQClimit.prop$red_zone),]
  plate.insideQClimit.prop <-
    reshape2::melt(plate.insideQClimit.prop, id.vars=c("experiment", "plate"))
  plate.insideQClimit.prop
}
