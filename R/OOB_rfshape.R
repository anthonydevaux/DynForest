#' OOB for random Forest
#'
#' @param rf
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
#' @import stringr
#' @import kmlShape
#' @import Evomorph
#' @import geomorph
#' @import pec
#' @import parallel
#' @import doParallel
#'
#' @keywords internal
OOB.rfshape <- function(rf, Curve=NULL, Scalar=NULL, Factor=NULL, Y, timeScale=0.1, d_out=0.1,
                        IBS.min = 0, IBS.max = NULL, cause = 1){

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

    allTimes <- sort(unique(c(0,Y$Y[,1])))

    if (is.null(IBS.max)){
      IBS.max <- max(allTimes)
    }

    Y.surv <- data.frame(time.event = Y$Y[order(Y$Y[,1]),1],
                         event = ifelse(Y$Y[order(Y$Y[,1]),2]==cause,1,0))

    # IPCW using all data
    ipcw.res <- pec::ipcw(formula = Surv(time.event, event) ~ 1,
                          data = Y.surv,
                          method = "marginal",
                          times=allTimes,
                          subjectTimes=allTimes)$IPCW.times

    comb <- function(x, ...) {
      mapply(rbind,x,...,SIMPLIFY=FALSE)
    }

    cl <- parallel::makeCluster(4)
    doParallel::registerDoParallel(cl)

    res.oob <- foreach(i=1:length(Y$id),
                       .combine='comb', .multicombine = TRUE,
                       .packages = c("pec", "prodlim")) %dopar%
      {
    # for (i in 1:length(Y$id)){
      indiv <- Y$id[i]
      w_y <- which(Y$id==indiv)
      #Y.surv = data.frame(time.event = Y$Y[w_y,1], event = Y$Y[w_y,2])
      pred.mat <- matrix(NA, nrow = ncol(rf$rf), ncol = length(allTimes))

      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(Y$id,BOOT)
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

          pred.node <- pred.MMT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale = timeScale)

          if (is.na(pred.node)){
            pred.mat[t,] <- NA
          }else{
            pred.mat[t,] <- rf$rf[,t]$Y_pred[[pred.node]]$traj
          }

        }
      }

      oob.pred <- apply(pred.mat, 2, mean, na.rm = TRUE)

      # individual IBS with IPCW using all data

      Di <- ifelse(Y$Y[w_y,1] <= allTimes, 1, 0) # D(t)
      pec.res <- list()
      pec.res$AppErr$matrix <- ipcw.res*(Di-oob.pred)^2 # BS(t)
      pec.res$models <- "matrix"
      pec.res$time <- allTimes
      pec.res$start <- 0
      pec.res$maxtime <- max(allTimes)
      class(pec.res) <- "pec"

      err <- pec::ibs(pec.res, start = IBS.min, times = IBS.max)[1] # IBS

      # pec.res <- pec::pec(object = t(oob.pred),
      #                     formula = Surv(time.event, event) ~ 1,
      #                     data = Y.surv, cens.model = "marginal",
      #                     exact = FALSE, times = allTimes,
      #                     maxtime = max(allTimes),
      #                     reference = FALSE)
      #
      # err <- pec::ibs(pec.res, start = 0, times = max(allTimes))[1]

      return(list(err=err,oob.pred=oob.pred))
    }

    parallel::stopCluster(cl)
    return(list(err=as.vector(res.oob$err),oob.pred=res.oob$oob.pred))
  }

  if (Y$type=="curve"){
    oob.pred <- list()
    #errdp <- rep(NA,length(unique(id)))

    for (i in 1:length(unique(Y$id))){
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
