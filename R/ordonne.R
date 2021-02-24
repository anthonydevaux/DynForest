#' Ordonne
#'
#' @param X
#' @param time
#' @param id
#'
#'
#' @keywords internal
ordonne <- function(X , time , id){
  mat  <- matrix(NA, length(unique(id)), length(unique(time)))
  for( i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    t_w <- time[w]
    w_time <- NULL
    for (j in 1:length(w)){
      w_time <- c(w_time, which(unique(time)==t_w[j]))
    }
    mat[i,w_time] <- X[w]
  }
  return(mat)
}
