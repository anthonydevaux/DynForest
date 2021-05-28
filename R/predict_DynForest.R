#' Predict with Dynamic random forests
#'
#' @param object : Dynamic random forest
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param timeScale [numeric]:
#' @param predTimes [numeric]:
#' @param ncores [numeric]:
#' @param ... : optional parameters to be passed to the low level function
#'
#' @import kmlShape
#' @import stringr
#' @import RiemBase
#' @import Evomorph
#' @import geomorph
#' @import parallel
#' @import doParallel
#'
#' @return
#' @export
#'
predict.DynForest <- function(object, Curve=NULL,Scalar=NULL,Factor=NULL, timeScale=0.1,
                              predTimes = NULL, t0 = NULL,
                              ncores = NULL, ...){

  ##########
  # Checking

  if (object$type=="surv"){

    if (is.null(t0)){
      stop("t0 value is needed for dynamic prediction !")
    }
    if (!is.null(predTimes)){
      if (all(predTimes<=t0)){
        stop("predTimes values should be greater than t0 time !")
      }
      if (any(predTimes<=t0)){
        warning("Only predTimes values greater than t0 time will be computed !")
        predTimes <- predTimes[!(predTimes <= t0)]
      }
    }else{
      predTimes <- object$times
    }

  }

  ##########

  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time,
                  model=object$Curve.model)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }

  if(is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  ## Puis on prend les prédicteurs:

  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs
  # On va les lires en mettant la maj sur les différents éléments qui le constituent :

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  Id.pred <- unique(get(Inputs[1])$id)
  pred.feuille <- matrix(0, ncol(object$rf), length(Id.pred))

  if (object$type=="factor"){
    pred.feuille <- as.data.frame(matrix(0, ncol(object$rf), length(Id.pred)))
  }

  # leaf predictions of new subjects
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  pred.feuille <- foreach(t=1:ncol(object$rf),
                       .combine='rbind', .multicombine = TRUE
                       #, .packages = c()
                       ) %dopar%
    {

      return(pred.MMT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor, timeScale))

    }

  parallel::stopCluster(cl)

  # for (t in 1:ncol(object$rf)){
  #   pred.feuille[t,] <- pred.MMT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor, timeScale)
  # }

  if (object$type=="scalar"){
    pred <- apply(pred.feuille, 2, "mean")
    return(pred)
  }

  if (object$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=object$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/ncol(object$rf)
    }
    prediction <- data.frame(pred=val, prob=proba)
    return(prediction)
  }

  if (object$type=="curve"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[,k]$Y_pred[[pred.feuille[k,l]]])[1]),object$rf[,k]$Y_pred[[pred.feuille[k,l]]]))
      }
      predouille <- kmlShape::meanFrechet(pred_courant, timeScale = timeScale)
      predouille <- cbind(predouille, rep(Id.pred[l],dim(predouille)[1]))
      pred <- rbind(pred, predouille)
    }
    names(pred) <- c("times", "traj", "ID")
  }

  if (object$type=="surv"){

    allTimes <- object$times

    # if (is.null(predTimes)){
    #   predTimes <- allTimes
    # }

    id.predTimes <- sapply(predTimes, function(x){ sum(allTimes <= x) })
    pred <- matrix(NA, nrow = length(Id.pred), ncol = length(predTimes))

    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- matrix(NA, nrow = ncol(object$rf), ncol = length(predTimes))
      for(k in 1:dim(pred.feuille)[1]){

        if (!is.na(pred.feuille[k,l])){
          pred_courant[k,] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]$traj[id.predTimes]
        }else{
          pred_courant[k,] <- NA
        }
      }

      pred[l,] <- apply(pred_courant, 2, mean, na.rm = TRUE)

    }
  }
  return(pred)
}
