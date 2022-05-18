#' Permute courbes among the subjects
#'
#' @param Courbes Numeric vector of the marker
#' @param id Numeric vector of id
#'
#' @keywords internal
permutation_courbes <- function(Courbes,id){
  perm <- sample(unique(id), length(unique(id))) #### on change les identifiants ::
  new <- NULL
  for (i in perm){
    w <- which(id==i)
    new <- c(new,Courbes[w])
  }
  return(new)
}
