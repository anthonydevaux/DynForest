#' Compute the impurity of a given vector
#'
#' @param Y Outcome data
#'
#' @importFrom DescTools Entropy
#'
#' @keywords internal
impurity <- function(Y){

  if (Y$type=="numeric"){
    if (length(Y$Y)==1){
      return(0)
    }
    return(var(Y$Y))
  }

  if (Y$type=="factor"){
    return(Entropy(table(Y$Y)))
  }

}
