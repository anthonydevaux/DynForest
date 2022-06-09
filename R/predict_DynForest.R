#' Prediction using dynamic random forests
#'
#' @param object \code{DynForest} object containing the dynamic random forest used on train data
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param predTimes (Only with survival outcome) Horizon times from landmark time \code{t0}
#' @param t0 Landmark time
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#' @param parallel Allow paralleling. Default value is TRUE.
#' @param ... optional parameters to be passed to the low level function
#'
#' @import stringr
#' @import parallel
#' @import doParallel
#'
#' @return Return the outcome of interest for the new subjects: matrix of probability of event of interest in survival mode, average value in regression mode and most likely value in classification mode
#'
#' @examples
#' \dontrun{
#' data(pbc2)
#'
#' # Build survival data
#' pbc2_surv <- unique(pbc2[,c("id","age","drug","sex","years","event")])
#'
#' # Define time-independent continuous covariate
#' cont_covar <- list(X = pbc2_surv[,"age", drop = FALSE],
#'                    id = pbc2_surv$id)
#'
#' # Define time-independent non continuous covariates
#' fact_covar <- list(X = pbc2_surv[,c("drug","sex")],
#'                    id = pbc2_surv$id)
#'
#' # Define time-dependent continuous markers
#' cont_traj <- list(X = pbc2[,c("serBilir","SGOT","albumin","alkaline")],
#'                   id = pbc2$id,
#'                   time = pbc2$time,
#'                   model = list(serBilir = list(fixed = serBilir ~ time,
#'                                                random = ~ time),
#'                                SGOT = list(fixed = SGOT ~ time + I(time^2),
#'                                            random = ~ time + I(time^2)),
#'                                albumin = list(fixed = albumin ~ time,
#'                                               random = ~ time),
#'                                alkaline = list(fixed = alkaline ~ time,
#'                                                random = ~ time))
#' )
#'
#' # Define outcome (survival here)
#' Y <- list(type = "surv",
#'           Y = Surv(pbc2_surv$years, factor(pbc2_surv$event)),
#'           id = pbc2_surv$id)
#'
#' # Run DynForest function
#' res_dyn <- DynForest(Curve = cont_traj, Factor = fact_covar, Scalar = cont_covar,
#'                      Y = Y, ntree = 200, imp = TRUE,
#'                      imp.group = list(group1 = c("serBilir","SGOT"),
#'                                       group2 = c("albumin","alkaline")),
#'                      mtry = 3, nodesize = 2, minsplit = 3,
#'                      cause = 2, seed = 1234)
#'
#' # Predict on new subjects using DynForest estimation (res_dyn object) from DynForest() function
#' pred_dyn <- predict(object = res_dyn,
#'                     Curve = cont_traj, Factor = fact_covar, Scalar = cont_covar,
#'                     predTimes = c(7,8), t0 = 4)
#'
#' }
#'
#' @export
predict.DynForest <- function(object, Curve = NULL, Scalar = NULL, Factor = NULL,
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
    Curve <- list(type="Curve",X=Curve$X,id=Curve$id,time=Curve$time,
                  model=object$Curve.model)
    if (!is.null(t0)){
      Curve$X <- Curve$X[which(Curve$time<=t0),]
      Curve$id <- Curve$id[which(Curve$time<=t0)]
      Curve$time <- Curve$time[which(Curve$time<=t0)]
    }
  }

  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="Scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="Factor",X=Factor$X,id=Factor$id)
  }

  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  #

  Inputs <- read.Xarg(c(Curve,Scalar,Factor))

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

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    pred.feuille <- foreach(t=1:ncol(object$rf),
                            .combine='rbind', .multicombine = TRUE
                            #, .packages = c()
    ) %dopar%
      {

        return(pred.MMT(object$rf[,t], Curve = Curve, Scalar = Scalar, Factor = Factor))

      }

    parallel::stopCluster(cl)

  }else{

    for (t in 1:ncol(object$rf)){
      #cat(t,"\n")
      pred.feuille[t,] <- pred.MMT(object$rf[,t], Curve = Curve, Scalar = Scalar,
                                   Factor = Factor)
    }

  }

  if (object$type=="scalar"){
    pred <- apply(pred.feuille, 2, "mean", na.rm = TRUE)
    return(pred)
  }

  if (object$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=object$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/sum(pred.all[[k]])
    }
    prediction <- data.frame(pred=val, prob=proba)
    return(prediction)
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
    pred.cause <- apply(pred[[as.character(object$cause)]][,-1, drop = FALSE],
                        MARGIN = 2,
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
