#' Title
#'
#' @param Shapes
#' @param id
#'
#' @keywords internal
permutation_shapes <- function(Shapes, id){
  perm <- sample(id,length(id))
  new <- array(0,dim=dim(Shapes)[1:3])
  for (i in 1:length(id)){
    new[,,i] <- Shapes[,,which(id==perm[i])]
  }
  return(new)
}
