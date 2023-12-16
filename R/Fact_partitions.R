#' Factor partitions finder
#'
#' This function is used to find all the unique partitions of k factors into 2 groups
#'
#' @param Factor Character vector
#' @param id List of id
#'
#' @keywords internal
Fact.partitions <- function(Factor, id){

  U <- unique(Factor)
  P <- Part.facts[[length(U)]]
  L <- vector("list", nrow(P))

  for (k in seq_len(nrow(P))){
    U_courant <- U[which(P[k,]==0)]
    W <- unlist(lapply(U_courant, function(m) which(Factor==m)))
    L[[k]] <- id[W]
  }
  return(L)
}
