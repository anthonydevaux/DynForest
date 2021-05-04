#' OOB tree
#'
#' @param tree
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Y
#' @param timeScale
#' @param d_out
#' @param IBS.min
#' @param IBS.max
#' @param cause
#'
#' @import kmlShape
#' @import Evomorph
#' @import stringr
#' @import RiemBase
#' @import pec
#' @import prodlim
#'
#' @keywords internal
OOB.tree <- function(tree, Curve=NULL, Scalar=NULL, Factor=NULL, Y, timeScale=0.1, d_out=0.1,
                     IBS.min = 0, IBS.max = NULL, cause = 1){

  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  BOOT <- tree$boot
  OOB <- setdiff(unique(Y$id), BOOT)
  xerror <- rep(NA,length(OOB))
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Curve_courant <- NULL

  if (Y$type=="curve" || Y$type=="surv"){
    if (Y$type=="surv"){
      allTimes <- sort(unique(c(0,Y$Y[,1])))

      if (is.null(IBS.max)){
        IBS.max <- max(allTimes)
      }

      Y.surv <- data.frame(time.event = Y$Y[order(Y$Y[,1]),1],
                           event = ifelse(Y$Y[order(Y$Y[,1]),2]==0,0,1))
      # IPCW using all data
      ipcw.res <- pec::ipcw(formula = Surv(time.event, event) ~ 1,
                            data = Y.surv,
                            method = "marginal",
                            times=allTimes,
                            subjectTimes=allTimes)$IPCW.times

    }

    for (i in OOB){

      id_wY <- which(Y$id== i)
      if (is.element("curve",inputs)==TRUE) {
        id_wXCurve <- which(Curve$id==i)
        Curve_courant <- list(type="curve",X=Curve$X[id_wXCurve,,drop=FALSE], id=Curve$id[id_wXCurve],time=Curve$time[id_wXCurve],
                              model=Curve$model)
      }

      if (is.element("factor",inputs)==TRUE){
        id_wXFactor <- which(Factor$id==i)
        Factor_courant <- list(type="factor",X=Factor$X[id_wXFactor,,drop=FALSE], id=Factor$id[id_wXFactor])
      }

      if (is.element("scalar",inputs)==TRUE){
        id_wXScalar <- which(Scalar$id==i)
        Scalar_courant <- list(type="scalar",X=Scalar$X[id_wXScalar,,drop=FALSE], id=Scalar$id[id_wXScalar])
      }

      pred_courant <- pred.MMT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale=timeScale)

      if (is.na(pred_courant)){

        xerror[which(OOB==i)] <- NA

      }else{

        if (Y$type == "curve"){
          xerror[which(OOB==i)] <- kmlShape::distFrechet(tree$Y_pred[[pred_courant]]$times, tree$Y_pred[[pred_courant]]$traj, Y$time[id_wY], Y$Y[id_wY], timeScale = d_out)^2
        }

        if (Y$type == "surv"){

          # individual IBS with IPCW using all data

          Di <- ifelse(Y$Y[id_wY,1] <= allTimes, 1, 0)*ifelse(Y$Y[id_wY,2]==cause,1,0) # D(t) = 1(s<Ti<s+t, event = cause)
          pec.res <- list()
          pec.res$AppErr$matrix <- ipcw.res*(Di-tree$Y_pred[[pred_courant]]$traj)^2 # BS(t)
          pec.res$models <- "matrix"
          pec.res$time <- allTimes
          pec.res$start <- 0
          pec.res$maxtime <- max(allTimes)
          class(pec.res) <- "pec"

          xerror[which(OOB==i)] <- pec::ibs(pec.res, start = IBS.min, times = IBS.max)[1] # IBS

        }

      }

    }

  }

  if (Y$type == "factor" || Y$type == "scalar"){
    w_XCurve <- NULL
    w_XScalar <- NULL
    w_XFactor <- NULL
    w_y <- NULL
    for (i in OOB){

      if (is.element("curve",inputs)==TRUE) w_XCurve <- c(w_XCurve, which(Curve$id==i))
      if (is.element("scalar",inputs)==TRUE) w_XScalar <- c(w_XScalar, which(Scalar$id==i))
      if (is.element("factor",inputs)==TRUE) w_XFactor <- c(w_XFactor, which(Factor$id==i))

      w_y <- c(w_y, which(Y$id==i))
    }

    if (is.element("curve",inputs)==TRUE) Curve_courant <- list(type="curve",X=Curve$X[w_XCurve,,drop=FALSE], id=Curve$id[w_XCurve],time=Curve$time[w_XCurve],
                                                                model=Curve$model)
    if (is.element("scalar",inputs)==TRUE) Scalar_courant  <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
    if (is.element("factor",inputs)==TRUE) Factor_courant  <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])

    pred <- pred.MMT(tree,Curve=Curve_courant,Scalar = Scalar_courant,Factor=Factor_courant)

    if (Y$type=="scalar"){xerror <- (Y$Y[w_y]-pred)^2}
    if (Y$type=="factor"){xerror <- 1*(pred!=Y$Y[w_y])}

  }

  return(mean(xerror, na.rm = T))
}
