#' Compute the Out-Of-Bag error (OOB error)
#'
#' @param DynForest_obj \code{DynForest} object containing the dynamic random forest used on train data
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#'
#' @return compute_OOBerror() function return a list with the following elements:\tabular{ll}{
#'    \code{data} \tab A list containing the data used to grow the trees \cr
#'    \tab \cr
#'    \code{rf} \tab A table with each tree in column. Provide multiple charactistics about the tree building \cr
#'    \tab \cr
#'    \code{type} \tab Outcome type \cr
#'    \tab \cr
#'    \code{times} \tab A numeric vector containing the time-to-event for all subjects \cr
#'    \tab \cr
#'    \code{cause} \tab Indicating the cause of interest \cr
#'    \tab \cr
#'    \code{causes} \tab A numeric vector containing the causes indicator \cr
#'    \tab \cr
#'    \code{Inputs} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Curve.model} \tab A list of longitudinal markers containing the formula used for modeling in the random forest \cr
#'    \tab \cr
#'    \code{param} \tab A list containing the hyperparameters \cr
#'    \tab \cr
#'    \code{xerror} \tab A numeric vector containing the OOB error for each tree \cr
#'    \tab \cr
#'    \code{oob.err} \tab A numeric vector containing the OOB error for each subject \cr
#'    \tab \cr
#'    \code{oob.pred} \tab Outcome prediction for all subjects \cr
#'    \tab \cr
#'    \code{IBS.range} \tab A vector containing the IBS min and max \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#' @export
#'
#' @examples
#' \dontrun{
#' # Compute OOB error
#' res_dyn_OOB <- compute_OOBerror(DynForest_obj = res_dyn)
#' }
compute_OOBerror <- function(DynForest_obj,
                             IBS.min = 0, IBS.max = NULL,
                             ncores = NULL){

  if (class(DynForest_obj)!="DynForest"){
    stop("'DynForest_obj' should be a 'DynForest' class!")
  }

  if (DynForest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(DynForest_obj$data$Y$Y[,1])
    }
  }

  rf <- DynForest_obj
  Curve <- rf$data$Curve
  Scalar <- rf$data$Scalar
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  ntree <- ncol(rf$rf)

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  ##############################

  print("OOB trees error...")
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  pck <- .packages()
  dir0 <- find.package()
  dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
  parallel::clusterExport(cl,list("pck","dir"),envir=environment())
  parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

  xerror <- pbsapply(1:ntree,
                     FUN=function(i){OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y,
                                              IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)},cl=cl)

  parallel::stopCluster(cl)

  # xerror <- rep(NA, ntree)
  # for (i in 1:ntree){
  #   cat(paste0("Tree ", i,"\n"))
  #   xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y,
  #                        IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)
  # }

  cat("OOB Forest error...")
  oob.err <- OOB.rfshape(rf, Curve = Curve, Scalar = Scalar, Factor = Factor, Y = Y,
                         IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause, ncores = ncores)
  cat("OK!\n")

  out <- list(data = rf$data,
              rf = rf$rf, type = rf$type, times = rf$times, cause = rf$cause, causes = rf$causes,
              Inputs = rf$Inputs, Curve.model = rf$Curve.model, param = rf$param,
              comput.time = rf$comput.time,
              xerror = xerror, oob.err = oob.err$err, oob.pred = oob.err$oob.pred,
              IBS.range = c(IBS.min, IBS.max))

  class(out) <- c("DynForest")

  return(out)

}
