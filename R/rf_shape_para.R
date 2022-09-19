#' Paralleled random survival Forest using multivariate longitudinal endogenous covariates
#'
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{curve}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error. Default is `NULL`.
#' @param ntree Number of trees to grow. Default value set to 200.
#' @param timeScale Experimental
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#' @param timeScale Experimental
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param seed Seed to replicate results
#' @param verbose A logical controlling the function progress. Default is \code{TRUE}
#'
#' @import foreach
#' @import doParallel
#' @import pbapply
#' @importFrom splines ns
#'
#' @keywords internal
rf_shape_para <- function(Curve = NULL, Scalar = NULL, Factor = NULL, Y, mtry, ntree, timeScale, ncores,
                          nsplit_option = "quantile", nodesize = 1, minsplit = 2, cause = 1,
                          seed = 1234, verbose = TRUE){

  if (!verbose){
    pbapply::pboptions(type="none")
  }else{
    pbapply::pboptions(type="timer")
  }

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  pck <- .packages()
  dir0 <- find.package()
  dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
  parallel::clusterExport(cl,list("pck","dir"),envir=environment())
  parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

  if (Y$type=="surv"){
    trees <- pbsapply(1:ntree, FUN=function(i){
      DynTree_surv(Y = Y, Curve=Curve,Scalar = Scalar,Factor = Factor, mtry = mtry,
                   nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause,
                   seed = seed*i)
    },cl=cl)
  }
  if (Y$type%in%c("factor","scalar")){
    trees <- pbsapply(1:ntree, FUN=function(i){
      DynTree(Y = Y, Curve=Curve,Scalar = Scalar,Factor = Factor, mtry = mtry,
                     nsplit_option = nsplit_option, nodesize = nodesize,
                     seed = seed*i)
    },cl=cl)
  }

  parallel::stopCluster(cl)

  # trees <- list()
  #
  # for (i in 1:ntree){
  #   cat(paste0("Tree ",i,"\n"))
  #
  #   #browser(expr = {i == 2})
  #
  #   if (Y$type=="surv"){
  #     trees[[i]] <- DynTree_surv(Y = Y, Curve = Curve, Scalar = Scalar, Factor = Factor,
  #                                mtry = mtry, nsplit_option = nsplit_option,
  #                                nodesize = nodesize, minsplit = minsplit, cause = cause,
  #                                seed = seed*i)
  #   }
  #   if (Y$type%in%c("factor","scalar")){
  #     trees[[i]] <- DynTree(Y = Y, Curve = Curve, Scalar = Scalar, Factor = Factor,
  #                           mtry = mtry, nsplit_option = nsplit_option,
  #                           nodesize = nodesize,
  #                           seed = seed*i)
  #   }
  # }

  return(trees)
}
