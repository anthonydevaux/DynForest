#' Compute the Out-Of-Bag error on the random survival forest
#'
#' @param rf Trees object resulting from \code{rf_shape_para} function
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{curve}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param timeScale Experimental
#' @param d_out Experimental
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
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
                        IBS.min = 0, IBS.max = NULL, cause = 1, ncores = NULL){

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

    Y.surv <- data.frame(id = Y$id,
                         time.event = Y$Y[,1],
                         event = ifelse(Y$Y[,2]==cause,1,0))

    Y.surv <- Y.surv[order(Y.surv$time.event, -Y.surv$event),]

    allTimes_IBS <- allTimes[which(allTimes>=IBS.min)]

    Y.surv <- Y.surv[which(Y.surv$time.event>=IBS.min),]

    # KM estimate of the survival function for the censoring
    G <- pec::ipcw(formula = Surv(time.event, event) ~ 1,
                   data = Y.surv,
                   method = "marginal",
                   times=allTimes_IBS,
                   subjectTimes=allTimes_IBS)

    OOB_IBS <- sort(Y.surv$id[which(Y.surv$time.event>=IBS.min)])

    comb <- function(x, ...) {
      mapply(rbind,x,...,SIMPLIFY=FALSE)
    }

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    res.oob <- foreach(i=1:length(OOB_IBS),
                       .combine='comb', .multicombine = TRUE
                       #,.packages = c("pec", "prodlim")
                       ) %dopar%
      {

    # for (i in 1:length(OOB_IBS)){

      indiv <- OOB_IBS[i]
      w_y <- which(Y$id==indiv)

      if (is.element("curve",inputs)==TRUE){
        w_XCurve <- which(Curve$id==indiv&Curve$time<=IBS.min) # only measurements until IBS.min
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

      pred.mat <- matrix(NA, nrow = ncol(rf$rf), ncol = length(allTimes_IBS))

      for (t in 1:ncol(rf$rf)){

        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(Y$id,BOOT)
        if (is.element(indiv, oob)== TRUE){

          pred.node <- tryCatch(pred.MMT(rf$rf[,t],Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale = timeScale),
                                         error = function(e) return(NA))

          if (is.na(pred.node)){
            pred.mat[t,] <- NA
          }else{
            pred.mat[t,] <- rf$rf[,t]$Y_pred[[pred.node]][[as.character(cause)]]$traj[which(allTimes>=IBS.min)]
          }

        }
      }

      oob.pred <- apply(pred.mat, 2, mean, na.rm = TRUE)

      # IPCW
      Wi_event <- (ifelse(Y$Y[w_y,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[w_y,2]!=0,1,0))/(G$IPCW.subjectTimes[which(Y.surv$id==indiv)])
      Wi_censored <- ifelse(Y$Y[w_y,1] > allTimes_IBS, 1, 0)/(G$IPCW.times)
      Wi <- Wi_event + Wi_censored

      # Individual Brier Score
      Di <- ifelse(Y$Y[w_y,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[w_y,2]==cause,1,0) # D(t) = 1(s<Ti<s+t, event = cause)
      pec.res <- list()
      pec.res$AppErr$matrix <- Wi*(Di-oob.pred)^2 # BS(t)
      pec.res$models <- "matrix"
      pec.res$time <- allTimes_IBS
      pec.res$start <- 0
      pec.res$maxtime <- max(allTimes_IBS)
      class(pec.res) <- "pec"

      #err <- pec::ibs(pec.res, start = IBS.min, times = IBS.max)[1] # IBS
      err <- pec::ibs(pec.res, start = IBS.min, times = max(allTimes))[1] # IBS

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
