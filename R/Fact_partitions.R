#' Factor partitions finder
#'
#' This function is used to find all the unique partitions of k factors into 2 groups
#'
#' @param Factor
#' @param id
#'
#' @keywords internal
Fact.partitions <- function(Factor, id){

  U <- unique(Factor)
  P <- Part.facts[[length(U)]]
  L <- list()
  for (k in 1:nrow(P)){
    w <- which(P[k,]==0)
    U_courant <- U[w]
    W <- NULL
    for (m in U_courant){
      W <- c(W,which(Factor==m))
    }
    L[[k]] <- id[W]
  }
  return(L)
}
