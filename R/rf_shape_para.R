#' Parallelized Dynamic random Forest
#'
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Shape
#' @param Image
#' @param Y
#' @param mtry
#' @param ntree
#' @param ncores
#' @param aligned.shape
#' @param timeScale
#' @param ...
#'
#' @import foreach
#' @import kmlShape
#' @import doParallel
#' @import pbapply
#' @importFrom splines ns
#'
#' @keywords internal
rf_shape_para <- function(Curve=NULL, Scalar=NULL, Factor=NULL,Shape=NULL,Image=NULL,Y,mtry,ntree, ncores, aligned.shape=FALSE,timeScale=0.1,
                          splitrule=NULL, nsplit_option=NULL, nodesize=1, ...){

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  trees <- pbsapply(1:ntree, FUN=function(i){
    Rtmax(Curve=Curve,Scalar = Scalar,Factor = Factor,Shape=Shape,Image=Image,Y,mtry, aligned.shape=aligned.shape,timeScale=timeScale,
          splitrule=splitrule, nsplit_option=nsplit_option, nodesize=nodesize, ...)
  },cl=cl)

  parallel::stopCluster(cl)

  return(trees)
}
