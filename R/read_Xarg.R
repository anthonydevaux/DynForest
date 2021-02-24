#' Read the parameters of the function
#'
#' @param z
#'
#'
#' @keywords internal
read.Xarg <- function(z){
  type <- NULL
  issou <- rep(NA, length(z))
  for (i in 1:length(z)){
    issou[i] <- is.null(z[i])
    if (issou[i]==FALSE) type <- c(type,z[i]$type)
  }
  return(type)
}
