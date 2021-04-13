#' Title
#'
#' @param pred
#' @param newtimes
#'
#' @return
#' @export
#'
#' @examples
combine_times <- function(pred, newtimes){

  newtimes.dt <- data.frame(times = newtimes, traj = NA)
  newpred <- rbind(pred, newtimes.dt)
  newpred <- newpred[order(newpred$times),]

  if (all(which(is.na(newpred$traj))!=1)){
    newpred$traj[which(is.na(newpred$traj))] <- newpred$traj[which(is.na(newpred$traj)) - 1]
  }else{
    newpred$traj[1] <- 1
    newpred$traj[which(is.na(newpred$traj))] <- newpred$traj[which(is.na(newpred$traj)) - 1]
  }

  return(newpred)

}
