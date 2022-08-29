#' Random forest with multivariate longitudinal endogenous covariates
#'
#' Built a random forest using multivariate longitudinal endogenous covariates
#'
#' @param timeData A data.frame containing the id and time measurements variables and the time-dependent predictors.
#' @param fixedData A data.frame containing the id variable and the time-fixed predictors. Non-continuous variables should be characterized as factor.
#' @param idVar A character indicating the name of variable to identify the subjects
#' @param timeVar A character indicating the name of time variable
#' @param timeVarModel A list for each time-dependent predictors containing a list of formula for fixed and random part from the mixed model
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{curve}", "\code{scalar}" or "\code{factor}"; .
#' @param ntree Number of trees to grow. Default value set to 200.
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error. Default is defined as the square root of the number of predictors.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#' @param seed Seed to replicate results
#' @param ... Optional parameters to be passed to the low level function
#'
#' @import stringr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import pbapply
#' @import survival
#' @import stats
#' @import utils
#'
#' @details The function currently supports survival (competing or single event), continuous or factor outcome.
#'
#' FUTUR IMPROVEMENTS:
#'
#' @return DynForest function return a list with the following elements:\tabular{ll}{
#'    \code{data} \tab A list containing the data used to grow the trees \cr
#'    \tab \cr
#'    \code{rf} \tab A table with each tree in column. Provide multiple charactistics about the tree building \cr
#'    \tab \cr
#'    \code{type} \tab Outcome type \cr
#'    \tab \cr
#'    \code{times} \tab A numeric vector containing the time-to-event for all subjects \cr
#'    \tab \cr
#'    \code{cause} \tab Indicating the cause of interest \cr
#'    \tab \cr
#'    \code{causes} \tab A numeric vector containing the causes indicator \cr
#'    \tab \cr
#'    \code{Inputs} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Curve.model} \tab A list of longitudinal markers containing the formula used for modeling in the random forest \cr
#'    \tab \cr
#'    \code{param} \tab A list containing the hyperparameters \cr
#'    \tab \cr
#'    \code{comput.time} \tab Computation time \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#'
#' @references Devaux et al., Random survival forests for competing risks with multivariate longitudinal endogenous covariates (2022), arXiv
#'
#' @seealso \code{\link{predict.DynForest} \link{plot_VIMP} \link{plot_gVIMP} \link{plot_mindepth}}
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
#' }
#' @export
#'
DynForest <- function(timeData = NULL, fixedData = NULL,
                      idVar = NULL, timeVar = NULL, timeVarModel = NULL,
                      Y = NULL, ntree = 200, mtry = NULL,
                      nodesize = 1, minsplit = 2, cause = 1,
                      OOB_error = TRUE, imp = FALSE, imp.group = NULL,
                      nsplit_option = "quantile",
                      IBS.min = 0, IBS.max = NULL, ncores = NULL,
                      seed = round(runif(1,0,10000)),
                      ...){

  debut <- Sys.time()

  if (is.null(mtry)){
    mtry <- round(sqrt(ifelse(!is.null(timeData), ncol(timeData)-2, 0) +
      ifelse(!is.null(fixedData), ncol(fixedData)-1, 0)))
  }

  # checking function
  checking(timeData = timeData, fixedData = fixedData,
           idVar = idVar, timeVar = timeVar, timeVarModel = timeVarModel, Y = Y,
           ntree = ntree, mtry = mtry, nodesize = nodesize, minsplit = minsplit)

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

  Inputs <- read.Xarg(c(Curve,Scalar,Factor))
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  # Output
  Y <- list(type = Y$type,
            id = Y$Y[,idVar],
            Y = Y$Y)

  if (Y$type=="surv"){
    Y$Y <- survival::Surv(Y$Y[,2], factor(Y$Y[,3]))
    Y$comp <- ifelse(length(unique(Y$Y[,2]))>2, TRUE, FALSE)
    causes <- sort(unique(Y$Y[which(Y$Y[,2]!=0),2]))
  }else{
    Y$Y <- subset(Y$Y, select = -get(idVar), drop = TRUE)
  }

  # gVIMP groups checking
  if (!is.null(imp.group)){
    if (!all(unlist(imp.group)%in%c(colnames(Curve$X),colnames(Factor$X),colnames(Scalar$X)))){
      stop("Unknown predictor(s) from 'imp.group' argument!")
    }
  }

  Importance <- gVIMP <- NULL
  oob.err <- xerror <- NULL

  # number of predictors
  nvar <- sum(sapply(Inputs, FUN = function(x) ncol(get(x)$X)))

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  ######### DynTree #########
  print("Growing Dynamic trees...")

  rf <-  rf_shape_para(Curve=Curve,Scalar=Scalar, Factor=Factor, Y=Y, mtry=mtry, ntree=ntree, ncores=ncores,
                       nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause, seed = seed)

  rf <- list(type=Y$type, rf=rf, levels=levels(Y$Y))

  ###########################

  cat("DynForest DONE!\n")

  if (Y$type == "surv"){
    drf <- list(data = list(Curve = Curve, Factor = Factor, Scalar = Scalar, Y = Y),
                rf = rf$rf, type = rf$type, times = sort(unique(c(0,Y$Y[,1]))), cause = cause, causes = causes,
                Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
                Curve.model = Curve$model, param = list(mtry = mtry, nodesize = nodesize,
                                                        minsplit = minsplit, ntree = ntree),
                comput.time = Sys.time() - debut)
    class(drf) <- c("DynForest")
    return(drf)
  }
  var.ini <- impurity(Y)
  varex <- 1 - mean(oob.err$err)/var.ini
  drf <- list(data = list(Curve = Curve, Factor = Factor, Scalar = Scalar, Y = Y),
              rf = rf$rf, type = rf$type, levels = rf$levels,
              varex = varex, Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
              Curve.model = Curve$model, param = list(mtry = mtry, nodesize = nodesize,
                                                      minsplit = NULL, ntree = ntree),
              comput.time = Sys.time() - debut)
  class(drf) <- c("DynForest")
  return(drf)
}
