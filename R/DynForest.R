#' Random forest with multivariate longitudinal endogenous covariates
#'
#' Build a random forest using multivariate longitudinal endogenous covariates
#'
#' @param timeData A data.frame containing the id and time measurements variables and the time-dependent predictors.
#' @param fixedData A data.frame containing the id variable and the time-fixed predictors. Categorical variables should be characterized as factor.
#' @param idVar A character indicating the name of variable to identify the subjects
#' @param timeVar A character indicating the name of time variable
#' @param timeVarModel A list for each time-dependent predictors containing a list of formula for fixed and random part from the mixed model
#' @param Y A list of output which should contain: \code{type} defines the nature of the outcome, can be "\code{surv}", "\code{numeric}" or "\code{factor}"; .
#' @param ntree Number of trees to grow. Default value set to 200.
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error. Default is defined as the square root of the number of predictors.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#' @param seed Seed to replicate results
#' @param verbose A logical controlling the function progress. Default is \code{TRUE}
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
#' @details The function currently supports survival (competing or single event), continuous or categorical outcome.
#'
#' FUTUR IMPLEMENTATIONS:
#' - Continuous longitudinal outcome
#'
#' @return DynForest function return a list with the following elements:\tabular{ll}{
#'    \code{data} \tab A list containing the data used to grow the trees \cr
#'    \tab \cr
#'    \code{rf} \tab A table with each tree in column. Provide multiple characteristics about the tree building \cr
#'    \tab \cr
#'    \code{type} \tab Outcome type \cr
#'    \tab \cr
#'    \code{times} \tab A numeric vector containing the time-to-event for all subjects \cr
#'    \tab \cr
#'    \code{cause} \tab Indicating the cause of interest \cr
#'    \tab \cr
#'    \code{causes} \tab A numeric vector containing the causes indicator \cr
#'    \tab \cr
#'    \code{Inputs} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Longitudinal.model} \tab A list of longitudinal markers containing the formula used for modeling in the random forest \cr
#'    \tab \cr
#'    \code{param} \tab A list containing the hyperparameters \cr
#'    \tab \cr
#'    \code{comput.time} \tab Computation time \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#'
#' @references Devaux A., Helmer C., Dufouil C., Genuer R., Proust-Lima C. (2022). Random survival forests for competing risks with multivariate longitudinal endogenous covariates. arXiv <doi: 10.48550/arXiv.2208.05801>
#'
#' @seealso \code{\link{summary.DynForest} \link{compute_OOBerror} \link{compute_VIMP} \link{compute_gVIMP} \link{predict.DynForest} \link{plot.DynForest}}
#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Get Gaussian distribution for longitudinal predictors
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
#'
#' # Sample 100 subjects
#' set.seed(1234)
#' id <- unique(pbc2$id)
#' id_sample <- sample(id, 100)
#' id_row <- which(pbc2$id%in%id_sample)
#'
#' pbc2_train <- pbc2[id_row,]
#'
#  Build longitudinal data
#' timeData_train <- pbc2_train[,c("id","time",
#'                                 "serBilir","SGOT",
#'                                 "albumin","alkaline")]
#'
#' # Create object with longitudinal association for each predictor
#' timeVarModel <- list(serBilir = list(fixed = serBilir ~ time,
#'                                      random = ~ time),
#'                      SGOT = list(fixed = SGOT ~ time + I(time^2),
#'                                  random = ~ time + I(time^2)),
#'                      albumin = list(fixed = albumin ~ time,
#'                                     random = ~ time),
#'                      alkaline = list(fixed = alkaline ~ time,
#'                                      random = ~ time))
#'
#' # Build fixed data
#' fixedData_train <- unique(pbc2_train[,c("id","age","drug","sex")])
#'
#' # Build outcome data
#' Y <- list(type = "surv",
#'           Y = unique(pbc2_train[,c("id","years","event")]))
#'
#' # Run DynForest function
#' res_dyn <- DynForest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#' }
#' @export
DynForest <- function(timeData = NULL, fixedData = NULL,
                      idVar = NULL, timeVar = NULL, timeVarModel = NULL,
                      Y = NULL, ntree = 200, mtry = NULL,
                      nodesize = 1, minsplit = 2, cause = 1,
                      nsplit_option = "quantile",
                      ncores = NULL,
                      seed = round(runif(1,0,10000)),
                      verbose = TRUE){

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
    Longitudinal <- list(type = "Longitudinal",
                         X = subset(timeData, select = -c(get(idVar), get(timeVar))),
                         id = timeData[,idVar],
                         time = timeData[,timeVar],
                         model = timeVarModel)
  }else{
    Longitudinal <- NULL
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
      Numeric <- list(type = "Numeric",
                      X = subset(fixedData, select = names(var_num[which(var_num==T)])),
                      id = fixedData[,idVar])
    }else{
      Numeric <- NULL
    }

  }else{
    Numeric <- NULL
    Factor <- NULL
  }

  Inputs <- read.Xarg(c(Longitudinal,Numeric,Factor))
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  # Outcome
  Y <- list(type = Y$type,
            id = Y$Y[,idVar],
            Y = Y$Y)

  if (Y$type=="surv"){

    causes <- sort(unique(Y$Y[which(Y$Y[,3]!=0),3]))

    if (length(unique(Y$Y[,3]))>2){
      Y$Y <- survival::Surv(Y$Y[,2], factor(Y$Y[,3]))
      Y$comp <- TRUE
    }else{
      Y$Y <- survival::Surv(Y$Y[,2], Y$Y[,3])
      Y$comp <- FALSE
    }
  }else{
    Y$Y <- subset(Y$Y, select = -get(idVar), drop = TRUE)
  }

  if (Y$type=="factor"){
    Ylevels <- unique(Y$Y)
  }else{
    Ylevels <- NULL
  }

  # number of predictors
  nvar <- sum(sapply(Inputs, FUN = function(x) ncol(get(x)$X)))

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  ######### DynTree #########

  rf <-  rf_shape_para(Longitudinal = Longitudinal, Numeric = Numeric, Factor = Factor,
                       timeVar = timeVar, Y = Y,
                       mtry = mtry, ntree = ntree, ncores = ncores,
                       nsplit_option = nsplit_option,
                       nodesize = nodesize, minsplit = minsplit,
                       cause = cause, seed = seed, verbose = verbose)

  rf <- list(type=Y$type, rf=rf)

  ###########################

  if (Y$type == "surv"){
    out <- list(data = list(Longitudinal = Longitudinal, Factor = Factor, Numeric = Numeric, Y = Y),
                rf = rf$rf, type = rf$type, timeVar = timeVar, times = sort(unique(c(0,Y$Y[,1]))), cause = cause, causes = causes,
                Inputs = list(Longitudinal = names(Longitudinal$X), Numeric = names(Numeric$X), Factor = names(Factor$X)),
                Longitudinal.model = Longitudinal$model, param = list(mtry = mtry, nodesize = nodesize,
                                                                      minsplit = minsplit, ntree = ntree),
                comput.time = Sys.time() - debut)
  }else{
    out <- list(data = list(Longitudinal = Longitudinal, Factor = Factor, Numeric = Numeric, Y = Y),
                rf = rf$rf, type = rf$type, timeVar = timeVar, levels = Ylevels,
                Inputs = list(Longitudinal = names(Longitudinal$X), Numeric = names(Numeric$X), Factor = names(Factor$X)),
                Longitudinal.model = Longitudinal$model, param = list(mtry = mtry, nodesize = nodesize,
                                                                      minsplit = NULL, ntree = ntree),
                comput.time = Sys.time() - debut)
  }

  class(out) <- c("DynForest")
  return(out)
}
