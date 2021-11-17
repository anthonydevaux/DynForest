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
                          nsplit_option = "quantile", nodesize = 1, minsplit = 2, cause = 1){

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,list("Y","minsplit","nodesize","Rtmax_surv"),
                          #envir = globalenv())
                          envir = environment())

  if (Y$type=="surv"){
    trees <- pbsapply(1:ntree, FUN=function(i){
      Rtmax_surv(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry,
                 nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause)
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
  #     #browser(expr = {i == 21})
  #     trees[[i]] <- Rtmax_surv(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry,
  #                              nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause)
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
