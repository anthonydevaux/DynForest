#' Prediction using dynamic random forests
#'
#' @param DynForest_obj \code{DynForest} or \code{DynForest_OOB} object containing the dynamic random forest used on train data
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
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
#'                     t0 = 4)
#'
#' }
#'
#' @export
predict.DynForest <- function(DynForest_obj,
                              timeData = NULL, fixedData = NULL,
                              idVar, timeVar,
                              t0 = NULL,
                              ncores = NULL, parallel = TRUE, ...){

  if (class(DynForest_obj)!="DynForest"){
    stop("'DynForest_obj' should be an object of 'DynForest' class!")
  }

  # checking function
  checking(DynForest_obj = DynForest_obj,
           timeData = timeData, fixedData = fixedData,
           idVar = idVar, timeVar = timeVar)

  # checking landmark/horizon times
  if (DynForest_obj$type=="surv"){

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

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
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
                  model = timeVarModel)
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
  pred.feuille <- matrix(0, ncol(DynForest_obj$rf), length(Id.pred))

  if (DynForest_obj$type=="factor"){
    pred.feuille <- as.data.frame(matrix(0, ncol(DynForest_obj$rf), length(Id.pred)))
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

    pred.feuille <- foreach(t=1:ncol(DynForest_obj$rf),
                            .combine='rbind', .multicombine = TRUE
                            #, .packages = c()
    ) %dopar%
      {

        return(pred.MMT(DynForest_obj$rf[,t], Curve = Curve, Scalar = Scalar, Factor = Factor))

      }

    parallel::stopCluster(cl)

  }else{

    for (t in 1:ncol(DynForest_obj$rf)){
      #cat(t,"\n")
      pred.feuille[t,] <- pred.MMT(DynForest_obj$rf[,t], Curve = Curve, Scalar = Scalar,
                                   Factor = Factor)
    }

  }

  if (DynForest_obj$type=="scalar"){
    pred_outcome <- apply(pred.feuille, 2, "mean", na.rm = TRUE)
    return(pred_outcome)
  }

  if (DynForest_obj$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=DynForest_obj$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/sum(pred.all[[k]])
    }
    pred_outcome <- data.frame(pred=val, prob=proba)
    return(pred_outcome)
  }

  if (DynForest_obj$type=="surv"){

    allTimes <- DynForest_obj$times
    predTimes <- c(t0, allTimes[which(allTimes>=t0)])

    id.predTimes <- sapply(predTimes, function(x){ sum(allTimes <= x) })

    pred <- lapply(DynForest_obj$causes, FUN = function(x){

      pred.cause <- matrix(NA, nrow = length(Id.pred), ncol = length(predTimes))
      rownames(pred.cause) <- Id.pred
      return(pred.cause)
    })
    names(pred) <- as.character(DynForest_obj$causes)

    for (l in 1:ncol(pred.feuille)){ # subject

      pred_courant <- lapply(DynForest_obj$causes, matrix, data = NA, nrow = ncol(DynForest_obj$rf), ncol = length(predTimes))
      names(pred_courant) <- as.character(DynForest_obj$causes)

      for (k in 1:nrow(pred.feuille)){ # tree

        if (!is.na(pred.feuille[k,l])){

          for (cause in as.character(DynForest_obj$causes)){
            if (!is.null(DynForest_obj$rf[,k]$Y_pred[[pred.feuille[k,l]]][[cause]]$traj[id.predTimes])){
              pred_courant[[cause]][k,] <- DynForest_obj$rf[,k]$Y_pred[[pred.feuille[k,l]]][[cause]]$traj[id.predTimes]
            }else{
              #pred_courant[[cause]][k,] <- rep(NA, length(id.predTimes))
              pred_courant[[cause]][k,] <- rep(0, length(id.predTimes))
            }
          }

        }else{
          for (cause in as.character(DynForest_obj$causes)){
            pred_courant[[cause]][k,] <- NA
          }
        }
      }

      for (cause in as.character(DynForest_obj$causes)){
        pred[[cause]][l,] <- apply(pred_courant[[cause]], 2, mean, na.rm = TRUE)
      }
    }

    # S landmark time / t horizon time
    # P(S<T<S+t|T>S) = ( P(T<S+t) - P(T<S) ) / P(T>S)
    #                = ( F(S+t) - F(S) ) / S(S)
    # A faire CR => S(S) n'est pas egale a 1-F(S) mais a la somme des Fj(S) avec j event
    pred_outcome <- apply(pred[[as.character(DynForest_obj$cause)]],
                        MARGIN = 2,
                        FUN = function(x) {
                          if (length(pred)>1){
                            surv <- 1 - Reduce("+", lapply(pred, FUN = function(x) x[,1]))
                          }else{
                            surv <- 1 - pred[[1]][,1]
                          }

                          return((x-pred[[as.character(DynForest_obj$cause)]][,1])/surv)
                        })

    output <- list(pred_outcome = pred_outcome,
                   pred_leaf = pred.feuille,
                   times = predTimes,
                   t0 = t0)

  }

  class(output) <- c("DynForestPred")
  return(output)

}
