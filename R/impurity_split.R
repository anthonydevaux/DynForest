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
impurity_split <- function(Y,split,cause=1){

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
        crr.res <- tryCatch(cmprsk::crr(ftime = Y$Y[,1], fstatus = Y$Y[,2], cov1 = split, failcode = cause),
                            error = function(e) return(list(converged = FALSE)))
        if (crr.res$converged){
          impur <- 2*pnorm(abs(crr.res$coef)/sqrt(diag(crr.res$var)), lower.tail=FALSE) # p-value (from emil package)
        }else{
          impur <- Inf
        }

        if (is.nan(impur)){
          impur <- Inf
        }

      }else{

        # logrank splitting rule
        vect_int <- c(3,6,12,24)
        nb_interval <- findInterval(sum(Y$Y[,2]), vect_int)+1
        idx_interval <- rep(1, length(Y$Y))
        random_interval <- 1
        browser()
        if (nb_interval > 1){
          print(paste0(nb_interval, " intervalles"))
          bornes <- quantile(Y$Y[,1][Y$Y[,2] == 1], probs = seq(0,1,length.out = nb_interval))
          idx_interval <- findInterval(Y$Y[,1], bornes)
          idx_interval <- ifelse(idx_interval == 0, 1, idx_interval) # pour changer 0 en 1
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
