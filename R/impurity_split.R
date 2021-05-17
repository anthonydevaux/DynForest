#' Impurity Split
#'
#' @param Y
#' @param split
#' @param cause
#'
#' @import survival
#' @importFrom cmprsk crr
#' @keywords internal
impurity_split <- function(Y,split,cause=1){

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
      if (Y$comp){

        # Fine & Gray splitting rule
        crr.res <- tryCatch(cmprsk::crr(ftime = Y$Y[,1], fstatus = Y$Y[,2], cov1 = split, failcode = cause),
                            error = function(e) return(list(converged = FALSE)))
        if (crr.res$converged){
          impur <- 2*pnorm(abs(crr.res$coef)/sqrt(diag(crr.res$var)), lower.tail=FALSE) # p-value (from emil package)
        }else{
          impur <- Inf
        }

      }else{

        # logrank splitting rule
        impur <- 1/(1+survival::survdiff(Y$Y~split)$chisq)

      }
      break
    }
  }
  return(list(impur=impur, imp_list=imp))
}
