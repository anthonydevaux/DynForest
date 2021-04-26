#' Parallelized Dynamic random Forest
#'
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Y
#' @param mtry
#' @param ntree
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
rf_shape_para <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Y, mtry, ntree, ncores, timeScale=0.1,
                          nsplit_option="quantile", nodesize=1, cause = 1){

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,list("Y"),envir=globalenv())

  trees <- pbsapply(1:ntree, FUN=function(i){
    Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry, timeScale = timeScale,
          nsplit_option = nsplit_option, nodesize = nodesize, cause = cause)
  },cl=cl)

  parallel::stopCluster(cl)

  # trees <- list()
  # for (i in 1:ntree){
  #   cat(paste0(i,"\n"))
  #   #browser(expr = {i == 46})
  #   trees[[i]] <- Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor, Y = Y, mtry = mtry, timeScale = timeScale,
  #                       nsplit_option = nsplit_option, nodesize = nodesize, cause = cause)
  # }

  return(trees)
}
