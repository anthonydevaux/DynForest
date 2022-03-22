#' Parallelized Dynamic random Forest
#'
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Y
#' @param mtry
#' @param ntree
#' @param timeScale
#' @param ncores
#' @param timeScale
#' @param nsplit_option
#' @param nodesize
#' @param cause
#' @param seed
#' @param ...
#'
#' @import foreach
#' @import kmlShape
#' @import doParallel
#' @import pbapply
#' @importFrom splines ns
#'
#' @keywords internal
rf_shape_para <- function(Curve = NULL, Scalar = NULL, Factor = NULL, Y, mtry, ntree, timeScale, ncores,
                          nsplit_option = "quantile", nodesize = 1, minsplit = 2, cause = 1,
                          seed = 1234){

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  pck <- .packages()
  dir0 <- find.package()
  dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
  parallel::clusterExport(cl,list("pck","dir"),envir=environment())
  parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

  if (Y$type=="surv"){
    trees <- pbsapply(1:ntree, FUN=function(i){
      Rtmax_surv(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry,
                 nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause,
                 seed = seed*i)
    },cl=cl)
  }else{
    trees <- pbsapply(1:ntree, FUN=function(i){
      Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry, timeScale = timeScale,
            nsplit_option = nsplit_option, nodesize = nodesize)
    },cl=cl)
  }

  parallel::stopCluster(cl)

  # trees <- list()
  #
  # if (Y$type=="surv"){
  #
  #   for (i in 1:ntree){
  #     cat(paste0("Tree ",i,"\n"))
  #     browser(expr = {i == 35})
  #     trees[[i]] <- Rtmax_surv(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry,
  #                              nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause, seed = seed*i)
  #   }
  #
  # }else{
  #
  #   for (i in 1:ntree){
  #     cat(paste0("Tree ",i,"\n"))
  #     #browser(expr = {i == 21})
  #     trees[[i]] <- Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry, timeScale = timeScale,
  #                         nsplit_option = nsplit_option, nodesize = nodesize)
  #   }
  #
  # }

  return(trees)
}
