#' Compute Out-Of-Bag error on the tree
#'
#' @param tree Tree object resulting from \code{Rtmax_surv} function
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{curve}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param timeScale Experimental
#' @param d_out Experimental
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param cause (Only with competing events) Number indicates the event of interest.
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

  Scalar_courant <- NULL
  Factor_courant <- NULL
  Curve_courant <- NULL

  if (Y$type=="curve" || Y$type=="surv"){

    if (Y$type=="surv"){ # survival outcome
      allTimes <- sort(unique(c(0,Y$Y[,1])))

      if (is.null(IBS.max)){
        IBS.max <- max(allTimes)
      }

      Y.surv <- data.frame(id = Y$id,
                           time.event = Y$Y[,1],
                           event = ifelse(Y$Y[,2]==cause,1,0))

      Y.surv <- Y.surv[order(Y.surv$time.event, -Y.surv$event),]

      #### NEW #####

      allTimes_IBS <- allTimes[which(allTimes>=IBS.min&allTimes<=IBS.max)]

      Y.surv <- Y.surv[which(Y.surv$time.event>=IBS.min),]

      # KM estimate of the survival function for the censoring
      G <- pec::ipcw(formula = Surv(time.event, event) ~ 1,
                     data = Y.surv,
                     method = "marginal",
                     times=allTimes_IBS,
                     subjectTimes=allTimes_IBS)

      OOB_IBS <- sort(Y.surv$id[which(Y.surv$id%in%OOB & Y.surv$time.event>=IBS.min)])

      xerror <- rep(NA,length(OOB_IBS))

      for (i in OOB_IBS){

        id_wY <- which(Y$id==i)
        if (is.element("curve",inputs)==TRUE) {

          if (IBS.min==0){
            id_wXCurve <- which(Curve$id==i) # all measurements
          }else{
            id_wXCurve <- which(Curve$id==i&Curve$time<=IBS.min) # only measurements until IBS.min
          }

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

        pred_courant <- tryCatch(pred.MMT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale=timeScale),
                                 error = function(e) return(NA)) # handle permutation issue

        if (is.na(pred_courant)){

          xerror[which(OOB_IBS==i)] <- NA

        }else{

          # IPCW
          Wi_event <- (ifelse(Y$Y[id_wY,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[id_wY,2]!=0,1,0))/(G$IPCW.subjectTimes[which(Y.surv$id==i)])
          Wi_censored <- ifelse(Y$Y[id_wY,1] > allTimes_IBS, 1, 0)/(G$IPCW.times)
          Wi <- Wi_event + Wi_censored

          # Individual Brier Score
          Di <- ifelse(Y$Y[id_wY,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[id_wY,2]==cause,1,0) # D(t) = 1(s<Ti<s+t, event = cause)
          pec.res <- list()
          pec.res$AppErr$matrix <- Wi*(Di-tree$Y_pred[[pred_courant]][[as.character(cause)]]$traj[which(allTimes%in%allTimes_IBS)])^2 # BS(t)
          pec.res$models <- "matrix"
          pec.res$time <- allTimes_IBS
          pec.res$start <- 0
          pec.res$maxtime <- max(allTimes_IBS)
          class(pec.res) <- "pec"

          #xerror[which(OOB==i)] <- pec::ibs(pec.res, start = IBS.min, times = IBS.max)[1] # IBS
          xerror[which(OOB_IBS==i)] <- pec::ibs(pec.res, start = IBS.min, times = max(allTimes_IBS))[1] # IBS

        }

      }

      ##############

    }


    if (Y$type=="curve"){ # curve outcome

      xerror <- rep(NA,length(OOB))

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

        pred_courant <- tryCatch(pred.MMT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant, timeScale=timeScale),
                                 error = function(e) return(NA)) # handle permutation issue

        if (is.na(pred_courant)){

          xerror[which(OOB==i)] <- NA

        }else{

          xerror[which(OOB==i)] <- kmlShape::distFrechet(tree$Y_pred[[pred_courant]]$times,
                                                         tree$Y_pred[[pred_courant]]$traj,
                                                         Y$time[id_wY],
                                                         Y$Y[id_wY],
                                                         timeScale = d_out)^2

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
