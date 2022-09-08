#' Prediction using dynamic random forests
#'
#' @param object \code{DynForest} or \code{DynForest_OOB} object containing the dynamic random forest used on train data
#' @param timeData A data.frame containing the id and time measurements variables and the time-dependent predictors.
#' @param fixedData A data.frame containing the id variable and the time-fixed predictors. Non-continuous variables should be characterized as factor.
#' @param idVar A character indicating the name of variable to identify the subjects
#' @param timeVar A character indicating the name of time variable
#' @param t0 Landmark time
#' @param ... Optional parameters to be passed to the low level function
#'
#' @import stringr
#' @importFrom methods is
#'
#' @return Return the outcome of interest for the new subjects: matrix of probability of event of interest in survival mode, average value in regression mode and most likely value in classification mode
#'
#' @examples
#' \dontrun{
#' # Predict on subjects still at risk at landmark time at 4 years
#' id_pred <- unique(pbc2_pred$id[which(pbc2_pred$years>4)])
#' pbc2_pred <- pbc2_pred[which(pbc2_pred$id%in%id_pred),]
#' timeData_pred <- pbc2_pred[,c("id", "time", "serBilir", "SGOT", "albumin", "alkaline")]
#' fixedData_pred <- unique(pbc2_pred[,c("id","age","drug","sex")])
#'
#' pred_dyn <- predict(object = res_dyn,
#'                     timeData = timeData_pred, fixedData = fixedData_pred,
#'                     idVar = "id", timeVar = "time",
#'                     t0 = 4)
#' }
#' @rdname predict.DynForest
#' @export
predict.DynForest <- function(object,
                              timeData = NULL, fixedData = NULL,
                              idVar, timeVar,
                              t0 = NULL, ...){

  if (!methods::is(object,"DynForest")){
    stop("'object' should be an object of 'DynForest' class!")
  }

  # checking function
  checking(DynForest_obj = object,
           timeData = timeData, fixedData = fixedData,
           idVar = idVar, timeVar = timeVar)

  # checking landmark/horizon times
  if (object$type=="surv"){

    if (is.null(t0)){
      stop("t0 value is needed for dynamic prediction !")
    }

  }

  # Select data before landmark time
  if (is.null(timeData)==FALSE){
    if (!is.null(t0)){
      timeData <- timeData[which(timeData[,timeVar]<=t0),]
    }
  }

  #####################
  # Handle missing data

  Inputs <- NULL

  if (!is.null(timeData)){

    timeData_id_noNA_list <- lapply(colnames(subset(timeData,
                                                    select = -c(get(idVar),get(timeVar)))),
                                    FUN = function(x){
                                      df <- timeData[,c(idVar,x)]
                                      colnames(df) <- c("ID","var")
                                      aggregate(var ~ ID, data = df, FUN = length)[,1]
                                    })

    timeData_id_noNA <- Reduce(intersect, timeData_id_noNA_list)

    if (length(timeData_id_noNA)>0){

      timeData <- timeData[which(timeData[,idVar]%in%timeData_id_noNA),]

    }else{
      stop("No subject has at least one measurement by marker !")
    }

    Inputs <- c(Inputs, "timeData")
  }

  if (!is.null(fixedData)){

    fixedData_na_row <- which(rowSums(is.na(fixedData))>0)

    if (length(fixedData_na_row)>0){
      fixedData <- fixedData[-fixedData_na_row,]
    }

    Inputs <- c(Inputs, "fixedData")
  }

  # all idnoNA
  idnoNA <- Reduce(intersect, lapply(Inputs, FUN = function(x) return(unique(get(x)[,idVar]))))

  # Keep id with noNA
  if (!is.null(timeData)){
    timeData <- timeData[which(timeData[,idVar]%in%idnoNA),]
  }

  if (!is.null(fixedData)){
    fixedData <- fixedData[which(fixedData[,idVar]%in%idnoNA),]
  }

  #####################

  # Inputs
  if (!is.null(timeData)){
    Curve <- list(type = "Curve",
                  X = subset(timeData, select = -c(get(idVar), get(timeVar))),
                  id = timeData[,idVar],
                  time = timeData[,timeVar],
                  model = object$Curve.model)
  }else{
    Curve <- NULL
  }

  if (!is.null(fixedData)){

    var_fact <- sapply(subset(fixedData, select = -get(idVar)),
                       FUN = function(x) inherits(x, c("character","factor")))

    var_num <- sapply(subset(fixedData, select = -get(idVar)),
                      FUN = function(x) inherits(x, c("numeric","integer")))

    if (length(var_fact[which(var_fact==T)])>0){
      Factor <- list(type = "Factor",
                     X = subset(fixedData, select = names(var_fact[which(var_fact==T)])),
                     id = fixedData[,idVar])
    }else{
      Factor <- NULL
    }

    if (length(var_num[which(var_num==T)])>0){
      Scalar <- list(type = "Scalar",
                     X = subset(fixedData, select = names(var_num[which(var_num==T)])),
                     id = fixedData[,idVar])
    }else{
      Scalar <- NULL
    }

  }

  #####################

  Id.pred <- as.integer(idnoNA)

  if (object$type=="surv"){

    allTimes <- object$times
    predTimes <- c(t0, allTimes[which(allTimes>=t0)])

    id.predTimes <- sapply(predTimes, function(x){ sum(allTimes <= x) })

    pred <- lapply(object$causes, FUN = function(x){

      lapply(Id.pred, FUN = function(x){
        pred_tree <- matrix(NA, nrow = ncol(object$rf), ncol = length(predTimes))
      })

    })

    names(pred) <- as.character(object$causes)

  }else{

    pred <- matrix(0, ncol(object$rf), length(Id.pred))

  }

  pred_leaf <- matrix(0, ncol(object$rf), length(Id.pred))

  ####################################
  # Leaf predictions of new subjects

  for (t in 1:ncol(object$rf)){

    pred_leaf[t,] <- pred.MMT(object$rf[,t],
                              Curve = Curve, Scalar = Scalar, Factor = Factor)

    if (object$type=="surv"){

      for (cause in as.character(object$causes)){

        for (indiv in seq(length(pred_leaf[t,]))){

          i.leaf <- pred_leaf[t,][indiv]

          pred_leaf_indiv <- object$rf[,t]$Y_pred[[i.leaf]][[cause]]$traj[id.predTimes]

          if (!is.null(pred_leaf_indiv)){
            pred[[cause]][[indiv]][t,] <- pred_leaf_indiv
          }else{
            pred[[cause]][[indiv]][t,] <- rep(0, length(id.predTimes))
          }

        }

      }

    }else{

      for (indiv in seq(length(pred_leaf[t,]))){

        i.leaf <- pred_leaf[t,indiv]

        pred_leaf_indiv <- object$rf[,t]$Y_pred[[i.leaf]]

        if (!is.null(pred_leaf_indiv)){
          pred[t,indiv] <- pred_leaf_indiv
        }else{
          pred[t,indiv] <- NA
        }

      }

    }

  }

  pred_out <- list(pred_leaf = pred_leaf,
                   pred = pred)

  if (object$type=="surv"){

    # Average CIF by subjects for each cause
    pred_cif_mean <- lapply(pred_out$pred, FUN = function(x){
      pred_cause_indiv <- t(sapply(x, FUN = function(y){
        apply(y, 2, mean, na.rm = TRUE)
      }))
      rownames(pred_cause_indiv) <- Id.pred
      return(pred_cause_indiv)
    })

    # S landmark time / t horizon time
    # P(S<T<S+t|T>S) = ( P(T<S+t) - P(T<S) ) / P(T>S)
    #                = ( F(S+t) - F(S) ) / S(S)
    # With competing risk S(S) = sum of Fj(S) avec j event
    pred_indiv <- apply(pred_cif_mean[[as.character(object$cause)]],
                        MARGIN = 2,
                        FUN = function(x) {
                          if (length(pred_cif_mean)>1){
                            surv <- 1 - Reduce("+", lapply(pred_cif_mean, FUN = function(x) x[,1]))
                          }else{
                            surv <- 1 - pred_cif_mean[[1]][,1]
                          }

                          return((x-pred_cif_mean[[as.character(object$cause)]][,1])/surv)
                        })

    output <- list(pred_indiv = pred_indiv,
                   pred_leaf = pred_out$pred_leaf,
                   times = predTimes,
                   t0 = t0)

  }

  if (object$type=="factor"){
    pred_indiv <- apply(pred_out$pred, 2, FUN = function(x) {
      tab_indiv <- table(x)
      return(names(which.max(tab_indiv)))
    })

    pred_indiv_proba <- apply(pred_out$pred, 2, FUN = function(x) {
      tab_indiv <- table(x)
      return(max(tab_indiv)/sum(tab_indiv))
    })

    names(pred_indiv) <- names(pred_indiv_proba) <- Id.pred

    output <- list(pred_indiv = pred_indiv,
                   pred_indiv_proba = pred_indiv_proba,
                   pred_indiv_tree = pred_out$pred,
                   pred_leaf = pred_out$pred_leaf,
                   t0 = t0)

  }

  if (object$type=="scalar"){
    pred_indiv <- apply(pred_out$pred, 2, "mean", na.rm = TRUE)

    names(pred_indiv) <- Id.pred

    output <- list(pred_indiv = pred_indiv,
                   pred_indiv_tree = pred_out$pred,
                   pred_leaf = pred_out$pred_leaf,
                   t0 = t0)
  }

  class(output) <- c("DynForestPred")
  return(output)

}
