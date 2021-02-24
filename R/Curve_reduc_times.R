#' Title
#'
#' @param time.init
#' @param traj.init
#' @param time.new
#'
#'
#' @keywords internal
Curve.reduc.times <- function(time.init , traj.init, time.new){
  new.Curve <- matrix(NA,length(time.new),2)
  for (j in 1:length(time.new)){
    w.time <- which.min(abs(time.new[j]-time.init))
    if (round(time.init[w.time]-time.new[j],5)==0){
      new.Curve[j,] <- c(time.new[j], traj.init[w.time])
    }
    else {
      t_g <- (time.new[j]>time.init[w.time])*(time.init[w.time]) + (time.new[j]<time.init[w.time])*(time.init[w.time-1])
      t_d <- (time.new[j]<time.init[w.time])*(time.init[w.time]) + (time.new[j]>time.init[w.time])*(time.init[w.time+1])
      Y_g <- (time.new[j]>time.init[w.time])*(traj.init[w.time]) + (time.new[j]<time.init[w.time])*(traj.init[w.time-1])
      Y_d <- (time.new[j]<time.init[w.time])*(traj.init[w.time]) + (time.new[j]>time.init[w.time])*(traj.init[w.time+1])
      d1 <- time.new[j]-t_g
      d2 <- t_d - time.new[j]
      new.Curve[j,] <- c(time.new[j], (1 - (d1/(d1+d2)))*Y_g + (1 - (d2/(d1+d2)))*Y_d)
    }
  }
  return(new.Curve)
}
