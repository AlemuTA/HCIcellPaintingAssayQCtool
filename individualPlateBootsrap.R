#' A function to estimate the multi-variate t-distribution for a particulat treatment (compound-concentration) and simulate new plates
#'
#' @param Y a data frame for summarized assay from the training plates in a wide format
#' @param AT a threshold, so that the mean of a feature value is set to be 0 if the aboslute value is less than AT.
#' By default, AT=1.
#' @param boot.size the number of bootsrap sampling (number of simulated plates). The default value is 2000.
#' @param settings.for.spcov a list of values for different values of the function \code{spcov()} in the 'spcov' package
#' @param seed an integer to control the reprodicibility of simulated plates
#'
#' @return a data frame
#' @import spcov spcov mnormt rmt
individualPlateBootsrap <- function(Y, AT=1, boot.size=2000,  seed=433,
                                    settings.for.spcov=list(step.size=10, eps=0.0001)){

  set.seed(seed)

  # calculate reguarized covariance matrix
  YY <- Y[,-c(1:3)]


  if(nrow(YY)>=ncol(YY)){
    Sig.hat <- list(Sigma=cov(YY))
  }else{
    # cross validation to obtaining the penality parameter (lambda)
    message("... ... running CV to obtain optimal regularizing parameter")
    landa.opt <- obtainOptimalLambda(Y=YY, settings.for.spcov = settings.for.spcov)
    #landa.opt <- 0.01 # for a quicker test

    # estimate covariance
    message("... ... obtaining regularized covariance estimate ")
    Sig.hat <- spcov(Sigma = diag(apply(YY, 2, var)),
                     S = var(YY)+diag(rep(settings.for.spcov$eps, ncol(YY))),
                     lambda = landa.opt,
                     step.size = settings.for.spcov$step.size,
                     trace = 0)
    #Sig.hat <- list(Sigma=diag(apply(YY, 2, var)))
  }

  message("... ... sampling new signatures")
  mean.vec <- apply(YY, 2, FUN = mean, na.rm=TRUE)
  mean.vec[abs(mean.vec)<AT] <- 0

  t.df = length(unique(Y$plate))
  YB <- mnormt::rmt(boot.size, mean = mean.vec,  S = Sig.hat$Sigma, df=t.df)
  YB <- as.data.frame(YB)
  colnames(YB) <- colnames(Y)[-c(1:3)]
  YB$plate <- paste0("boot_plate", seq_len(boot.size))
  YB$experiment <- paste0("boot_experiment", seq_len(boot.size))
  YB$compound_cnctr <- unique(Y$compound_cnctr)
  YB <- reshape2::melt(YB, id.vars=c("plate", "experiment", "compound_cnctr"))
  colnames(YB)[colnames(YB)=="variable"] <- "feature"
  YB
}
