#' Impurity Split
#'
#' @param Y
#' @param split
#'
#' @import survival
#' @keywords internal
impurity_split <- function(Y,split){
  impur <- 0
  imp <- list()
  for (i in 1:2){
    fils <- unique(Y$id)[which(split==i)]
    prop <- length(fils)/length(unique(Y$id))
    if (Y$type=="curve"){
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type="curve",Y=Y$Y[w],id=Y$id[w],time=Y$time[w]))
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type=="scalar" || Y$type=="factor") {
      w <- NULL
      for (j in 1:length(fils)){
        w <- c(w, which(Y$id==fils[j]))
      }
      imp[[i]] <- impurity(list(type=Y$type,Y=Y$Y[w],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }

    if (Y$type == "surv"){
      impur= 1/(1+survdiff(Y$Y~split)$chisq)
      break
    }
  }
  return(list(impur=impur, imp_list=imp))
}
