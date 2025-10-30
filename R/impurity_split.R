#' Impurity Split
#'
#' @param Y Outcome data
#' @param split Vector containing the subjects groups
#' @param cause (Only with competing events) Number indicates the event of interest.
#'
#' @import survival
#' @importFrom cmprsk crr
#'
#' @keywords internal
impurity_split <- function(Y,split,cause=1, randsplit = FALSE, splits_evt = NULL) {

  impur <- 0
  imp <- list()
  for (i in 1:2){
    fils <- unique(Y$id)[which(split==i)]
    prop <- length(fils)/length(unique(Y$id))

    if (Y$type=="numeric" || Y$type=="factor") {
      w <- which(Y$id%in%fils)
      imp[[i]] <- impurity(list(type=Y$type,Y=Y$Y[w],id=Y$id[w]))
      impur <- impur + imp[[i]]*prop
    }
    if (Y$type == "surv"){
      if (Y$comp){
        # Fine & Gray splitting rule
        if (randsplit){
          if (is.null(splits_evt)){
            splits_evt <- c(3,6,12,24)
          }
        } else {
          splits_evt <- Inf
        }
        nb_interval <- findInterval(sum(Y$Y[,2]==cause), splits_evt)+1
        idx_interval <- rep(1, length(Y$Y))
        random_interval <- 1
        if (nb_interval > 1){
          bornes <- quantile(Y$Y[,1][Y$Y[,2] == cause], probs = seq(0,1,length.out = nb_interval+1))
          idx_interval <- findInterval(Y$Y[,1], bornes)
          idx_interval <- ifelse(idx_interval == 0, 1, idx_interval)
          idx_interval <- ifelse(idx_interval == nb_interval+1, nb_interval, idx_interval)
          random_interval <- ceiling(runif(1,0,nb_interval))
        }

        crr.res <- tryCatch(cmprsk::crr(ftime = Y$Y[,1][idx_interval == random_interval], fstatus = Y$Y[,2][idx_interval == random_interval], cov1 = split[idx_interval == random_interval], failcode = cause),
                            error = function(e) return(list(converged = FALSE)))

        if (crr.res$converged & (length(Y$Y[,1]) == length(split))){ # condition on length to avoid troubles
          impur <- 2*pnorm(abs(crr.res$coef)/sqrt(diag(crr.res$var)), lower.tail=FALSE) # p-value (from emil package)
        }else{
          impur <- Inf
        }


        if (is.nan(impur)){
          impur <- Inf
        }

      } else {
        # logrank splitting rule
        if (randsplit){
          if (is.null(splits_evt)){
            splits_evt <- c(3,6,12,24)
          }
        } else {
          splits_evt <- Inf
        }

        nb_interval <- findInterval(sum(Y$Y[,2]), splits_evt)+1
        idx_interval <- rep(1, length(Y$Y))
        random_interval <- 1
        if (nb_interval > 1){
          bornes <- quantile(Y$Y[,1][Y$Y[,2] == 1], probs = seq(0,1,length.out = nb_interval+1))
          idx_interval <- findInterval(Y$Y[,1], bornes)
          idx_interval <- ifelse(idx_interval == 0, 1, idx_interval) # pour changer 0 en 1
          idx_interval <- ifelse(idx_interval == nb_interval+1, nb_interval, idx_interval)
          random_interval <- ceiling(runif(1,0,nb_interval))
        }
        surv.res <- tryCatch(survival::survdiff(Y$Y[idx_interval == random_interval]~split[idx_interval == random_interval]),
                             error = function(e) return(list(chisq = NULL)))

        if (!is.null(surv.res$chisq)){
          impur <- 1/(1+surv.res$chisq)
        }else{
          impur <- Inf
        }

      }
      break
    }
  }
  return(list(impur=impur, imp_list=imp))
}
