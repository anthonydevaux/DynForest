#' OOB for random Forest
#'
#' @param rf
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Y
#' @param timeScale
#' @param d_out
#'
#' @import stringr
#' @import kmlShape
#' @import Evomorph
#' @import geomorph
#'
#' @keywords internal
OOB.rfshape <- function(rf, Curve=NULL, Scalar=NULL, Factor=NULL, Y, timeScale=0.1, d_out=0.1){

  ### Pour optimiser le code il faudra virer cette ligne et ne le calculer qu'une seule fois !
  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs


  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }


  err <- rep(NA,length(unique(Y$id)))

  Curve_courant <- NULL
  Scalar_courant <- NULL
  Factor_courant <- NULL

  if (Y$type=="surv"){
    oob.pred <- list()
    #errdp <- rep(NA,length(unique(id)))

    for (i in 1:length(unique(Y$id))){
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      courbe2 = Y$Y[w_y]
      pred_courant <- rep(0, length(w_y))
      n_ind = 0
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){
          n_ind= n_ind+1

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve],
                                  model=Curve$model)
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          pred <- pred.MMT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale = timeScale)

          courbe <- rf$rf[,t]$Y_pred[[pred]] ## Il faut les mettre aux mêmes temps que Y$time[w_y]

          for (j in 1:length(w_y)){
            courbe2[j] = courbe[which.min(abs(Y$time[w_y][j]-courbe[,1])),2]
          }
          pred_courant <- pred_courant + courbe2
        }
      }
      oob.pred[[i]] <-  data.frame(times=Y$time[w_y], traj=pred_courant/n_ind)
      err[i] <- mean((Y$Y[w_y]-pred_courant/n_ind)^2)
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="curve"){
    oob.pred <- list()
    #errdp <- rep(NA,length(unique(id)))

    for (i in 1:length(unique(Y$id))){
      print(i)
      indiv <- unique(Y$id)[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve],
                                  model=Curve$model)
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          pred <- pred.MMT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale = timeScale)

          courbe <- rf$rf[,t]$Y_pred[[pred]]
          pred_courant <- rbind(cbind(rep(t,dim(courbe)[1]),courbe),pred_courant)
        }
      }
      mean_pred <- meanFrechet(pred_courant, timeScale = d_out)
      dp <- as.data.frame(Curve.reduc.times(mean_pred$times, mean_pred$traj, Y$time[w_y]))
      names(dp) <- c("x","y")
      oob.pred[[i]] <- dp
      err[i] <- distFrechet(dp$x, dp$y, Y$time[w_y], Y$Y[w_y], timeScale = d_out)^2
    }
    return(list(err=err,oob.pred=oob.pred))
  }

  if (Y$type=="scalar"){
    oob.pred <- rep(NA, length(unique(Y$id)))
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- NULL
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve],
                                  model=Curve$model)
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          pred <- pred.MMT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale = timeScale)
          pred_courant <- c(pred_courant, pred)
        }
      }
      oob.pred[i] <- mean(pred_courant)
      err[i] <- (oob.pred[i]-Y$Y[w_y])^2
    }
  }

  if (Y$type=="factor"){
    oob.pred <- factor(rep(NA, length(unique(Y$id))), levels=rf$levels)
    #errdp <- rep(NA,length(unique(id)))
    for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      pred_courant <- factor(rep(NA, length(unique(Y$id))), levels=rf$levels)
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id),BOOT)
        if (is.element(indiv, oob)== TRUE){

          if (is.element("curve",inputs)==TRUE){
            w_XCurve <- which(Curve$id== indiv)
            Curve_courant <- list(type="curve", X=Curve$X[w_XCurve,, drop=FALSE], id=Curve$id[w_XCurve], time=Curve$time[w_XCurve],
                                  model=Curve$model)
          }

          if (is.element("scalar",inputs)==TRUE){
            w_XScalar <- which(Scalar$id== indiv)
            Scalar_courant <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
          }

          if (is.element("factor",inputs)==TRUE){
            w_XFactor <- which(Factor$id== indiv)
            Factor_courant <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
          }

          pred <- pred.MMT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale = timeScale)
          pred_courant[t] <- pred
        }
      }
      pred_courant <- na.omit(pred_courant)
      oob.pred[i] <- as.factor(attributes(which.max(table(pred_courant))))
    }
    err <- 1*(oob.pred!=Y$Y)
  }

  return(list(err=err,oob.pred=oob.pred))
}
