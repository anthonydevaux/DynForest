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
                              ncores = NULL, parallel = TRUE, ...){

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
      predTimes <- object$times[!(object$times <= t0)]
    }

  }

  #############
  # t0 landmark

  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time,
                  model=object$Curve.model)
    if (!is.null(t0)){
      Curve$X <- Curve$X[which(Curve$time<=t0),]
      Curve$id <- Curve$id[which(Curve$time<=t0)]
      Curve$time <- Curve$time[which(Curve$time<=t0)]
    }
  }

  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }

  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  #

  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  #####################
  # Handle missing data

  if (!is.null(Curve)){ # Curve

    Curve_df <- cbind(ID = Curve$id, time = Curve$time, Curve$X)

    curve_id_noNA_list <- lapply(colnames(Curve_df)[3:ncol(Curve_df)], FUN = function(x){
      df <- Curve_df[,c("ID",x)]
      colnames(df) <- c("ID","var")
      aggregate(var ~ ID, data = df, FUN = length)[,1]
    })

    curve_id_noNA <- Reduce(intersect, curve_id_noNA_list)

    if (length(curve_id_noNA)>0){

      wCurve <- which(Curve$id%in%curve_id_noNA)

      Curve$X <- Curve$X[wCurve,, drop = FALSE]
      Curve$id <- Curve$id[wCurve]
      Curve$time <- Curve$time[wCurve]
    }else{
      stop("No subject has at least one measurement by marker !")
    }
  }

  if (!is.null(Scalar)){ # Scalar

    scalar_na_row <- which(rowSums(is.na(Scalar$X))>0)

    if (length(scalar_na_row)>0){
      Scalar$X <- Scalar$X[-scalar_na_row,, drop = FALSE]
      Scalar$id <- Scalar$id[-scalar_na_row]
    }
  }

  if (!is.null(Factor)){ # Factor

    factor_na_row <- which(rowSums(is.na(Factor$X))>0)

    if (length(factor_na_row)>0){
      Factor$X <- Factor$X[-factor_na_row,, drop = FALSE]
      Factor$id <- Factor$id[-factor_na_row]
    }
  }

  # all idnoNA
  idnoNA <- Reduce(intersect, lapply(Inputs, FUN = function(x) return(unique(get(x)$id))))

  # Keep id with noNA
  if (!is.null(Curve)){
    Curve$X <- Curve$X[which(Curve$id%in%idnoNA),, drop = FALSE]
    Curve$id <- Curve$id[which(Curve$id%in%idnoNA)]
    Curve$time <- Curve$time[which(Curve$id%in%idnoNA)]
  }

  if (!is.null(Scalar)){
    Scalar$X <- Scalar$X[which(Scalar$id%in%idnoNA),, drop = FALSE]
    Scalar$id <- Scalar$id[which(Scalar$id%in%idnoNA)]
  }

  if (!is.null(Factor)){
    Factor$X <- Factor$X[which(Factor$id%in%idnoNA),, drop = FALSE]
    Factor$id <- Factor$id[which(Factor$id%in%idnoNA)]
  }

  #####################

  Id.pred <- unique(get(Inputs[1])$id)
  pred.feuille <- matrix(0, ncol(object$rf), length(Id.pred))

  if (object$type=="factor"){
    pred.feuille <- as.data.frame(matrix(0, ncol(object$rf), length(Id.pred)))
  }

  # leaf predictions of new subjects

  if (parallel){

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

  }else{

    for (t in 1:ncol(object$rf)){
      #cat(t,"\n")
      pred.feuille[t,] <- pred.MMT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor, timeScale)
    }

  }

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
    for (l in 1:dim(pred.feuille)[2]){ # subject
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){ # tree
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
    predTimes <- c(t0, predTimes)

    id.predTimes <- sapply(predTimes, function(x){ sum(allTimes <= x) })

    pred <- lapply(object$causes, FUN = function(x){

      pred.cause <- matrix(NA, nrow = length(Id.pred), ncol = length(predTimes))
      rownames(pred.cause) <- Id.pred
      return(pred.cause)
    })
    names(pred) <- as.character(object$causes)

    for (l in 1:ncol(pred.feuille)){ # subject

      pred_courant <- lapply(object$causes, matrix, data = NA, nrow = ncol(object$rf), ncol = length(predTimes))
      names(pred_courant) <- as.character(object$causes)

      for(k in 1:nrow(pred.feuille)){ # tree

        if (!is.na(pred.feuille[k,l])){

          for (cause in as.character(object$causes)){
            if (!is.null(object$rf[,k]$Y_pred[[pred.feuille[k,l]]][[cause]]$traj[id.predTimes])){
              pred_courant[[cause]][k,] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]][[cause]]$traj[id.predTimes]
            }else{
              #pred_courant[[cause]][k,] <- rep(NA, length(id.predTimes))
              pred_courant[[cause]][k,] <- rep(0, length(id.predTimes))
            }
          }

        }else{
          for (cause in as.character(object$causes)){
            pred_courant[[cause]][k,] <- NA
          }
        }
      }

      for (cause in as.character(object$causes)){
        pred[[cause]][l,] <- apply(pred_courant[[cause]], 2, mean, na.rm = TRUE)
      }
    }

    # S landmark time / t horizon time
    # P(S<T<S+t|T>S) = ( P(T<S+t) - P(T<S) ) / P(T>S)
    #                = ( F(S+t) - F(S) ) / S(S)
    # A faire CR => S(S) n'est pas egale a 1-F(S) mais a la somme des Fj(S) avec j event
    pred.cause <- apply(pred[[as.character(object$cause)]][,-1], MARGIN = 2,
                        FUN = function(x) {
                          if (length(pred)>1){
                            surv <- 1 - Reduce("+", lapply(pred, FUN = function(x) x[,1]))
                          }else{
                            surv <- 1 - pred[[1]][,1]
                          }

                          return((x-pred[[as.character(object$cause)]][,1])/surv)
                        })

  }
  return(pred.cause)
}
