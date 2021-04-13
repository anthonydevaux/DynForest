#' Impurity
#'
#' Compute the impurity of a given vector
#'
#' @param Y
#' @param timeScale
#'
#' @import kmlShape
#' @import RiemBase
#' @import DescTools
#' @import Evomorph
#' @import geomorph
#' @import stats
#'
#'
#' @keywords internal
impurity <- function(Y, timeScale=0.1){

  if (Y$type=="curve"){
    traj <- Y$Y
    id <- Y$id
    time <- Y$time
    imp <- 0
    trajLong <- data.frame(id=id,time=time,traj=traj)
    meanF <- meanFrechet(trajLong = trajLong, timeScale = timeScale)
    for (i in unique(id)){
      imp <- imp + distFrechet(meanF$times, meanF$traj, time[which(id==i)], traj[which(id==i)], timeScale = timeScale)^2
    }
    imp <- imp/length(unique(id))
    return(imp)
  }

  if (Y$type=="scalar"){
    if (length(Y$Y)==1){
      return(0)
    }
    return(var(Y$Y))
  }

  if (Y$type=="factor"){
    return(Entropy(table(Y$Y)))
  }

}
