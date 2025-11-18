#' Extend predictions for new times
#'
#' @param pred Prediction matrix
#' @param newtimes Numeric vector of new times
#' @param type Character indicating survival function (\code{type}="surv") or risk function (\code{type}="risk")
#'
#' @importFrom zoo na.locf
#'
#' @keywords internal
#' @noRd
combine_times <- function(pred, newtimes, type = "surv"){

  newtimes <- unique(sort(c(0,newtimes, pred$times)))

  newtimes.dt <- data.frame(times = newtimes)

  newpred <- merge(newtimes.dt, pred, all.x = T)

  newpred$traj <- zoo::na.locf(newpred$traj, na.rm = F)

  if (any(is.na(newpred$traj))){
    if (type == "surv"){
      newpred$traj[which(is.na(newpred$traj))] <- 1
    }
    if (type == "risk"){
      newpred$traj[which(is.na(newpred$traj))] <- 0
    }
  }

  return(newpred)

}
