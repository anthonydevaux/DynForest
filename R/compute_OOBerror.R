#' Compute the Out-Of-Bag error (OOB error)
#'
#' @param DynForest_obj \code{DynForest} object containing the dynamic random forest used on train data
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#'
#' @importFrom methods is
#'
#' @return \code{compute_OOBerror()} function return a list with the following elements:\tabular{ll}{
#'    \code{data} \tab A list containing the data used to grow the trees \cr
#'    \tab \cr
#'    \code{rf} \tab A table with each tree in column. Provide multiple characteristics about the tree building \cr
#'    \tab \cr
#'    \code{type} \tab Outcome type \cr
#'    \tab \cr
#'    \code{times} \tab A numeric vector containing the time-to-event for all subjects \cr
#'    \tab \cr
#'    \code{cause} \tab Indicating the cause of interest \cr
#'    \tab \cr
#'    \code{causes} \tab A numeric vector containing the causes indicator \cr
#'    \tab \cr
#'    \code{Inputs} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Longitudinal.model} \tab A list of longitudinal markers containing the formula used for modeling in the random forest \cr
#'    \tab \cr
#'    \code{param} \tab A list containing the hyperparameters \cr
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
#' \donttest{
#' data(pbc2)
#'
#' # Get Gaussian distribution for longitudinal predictors
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
#'
#' # Sample 100 subjects
#' set.seed(1234)
#' id <- unique(pbc2$id)
#' id_sample <- sample(id, 100)
#' id_row <- which(pbc2$id%in%id_sample)
#'
#' pbc2_train <- pbc2[id_row,]
#'
#  Build longitudinal data
#' timeData_train <- pbc2_train[,c("id","time",
#'                                 "serBilir","SGOT",
#'                                 "albumin","alkaline")]
#'
#' # Create object with longitudinal association for each predictor
#' timeVarModel <- list(serBilir = list(fixed = serBilir ~ time,
#'                                      random = ~ time),
#'                      SGOT = list(fixed = SGOT ~ time + I(time^2),
#'                                  random = ~ time + I(time^2)),
#'                      albumin = list(fixed = albumin ~ time,
#'                                     random = ~ time),
#'                      alkaline = list(fixed = alkaline ~ time,
#'                                      random = ~ time))
#'
#' # Build fixed data
#' fixedData_train <- unique(pbc2_train[,c("id","age","drug","sex")])
#'
#' # Build outcome data
#' Y <- list(type = "surv",
#'           Y = unique(pbc2_train[,c("id","years","event")]))
#'
#' # Run DynForest function
#' res_dyn <- DynForest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Compute OOB error
#' res_dyn_OOB <- compute_OOBerror(DynForest_obj = res_dyn, ncores = 2)
#' }
compute_OOBerror <- function(DynForest_obj,
                             IBS.min = 0, IBS.max = NULL,
                             ncores = NULL){

  if (!methods::is(DynForest_obj,"DynForest")){
    stop("'DynForest_obj' should be a 'DynForest' class!")
  }

  if (DynForest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(DynForest_obj$data$Y$Y[,1])
    }
  }

  rf <- DynForest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  timeVar <- rf$timeVar
  Y <- rf$data$Y
  ntree <- ncol(rf$rf)

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  # Internal function to compute the OOB error for each subject
  oob.err <- OOB.rfshape(rf, Longitudinal = Longitudinal, Numeric = Numeric, Factor = Factor, Y = Y,
                         timeVar = timeVar, IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause,
                         ncores = ncores)

  out <- DynForest_obj
  out$oob.err <- oob.err$err
  out$oob.pred <- oob.err$oob.pred
  out$IBS.range <- c(IBS.min, IBS.max)

  class(out) <- c("DynForest")

  return(out)

}
