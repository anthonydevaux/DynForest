#' Predict Mixed Model Tree
#'
#' @param tree : Mixed Model Tree
#' @param Curve [list]: A list that contains the input curves.
#' @param Scalar [list]: A list that contains the input scalars.
#' @param Factor [list]: A list that contains the input factors.
#' @param timeScale [numeric]: Time scale for the input and output curves (\code{timeScale=0.1} by default)
#'
#' @import stringr
#' @import geomorph
#' @import kmlShape
#' @import Evomorph
#' @import RiemBase
#'
#' @return
#' @export
#'
pred.MMT <- function(tree, Curve=NULL,Scalar=NULL,Factor=NULL,timeScale=0.1){

  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  id.pred <- unique(get(Inputs[1])$id)

  if (tree$Y$type=="factor"){
    pred <- factor(rep(NA, length(id.pred)),levels=tree$Ylevels)
  }else{pred <- rep(NA,length(id.pred))}

  for (i in 1:length(id.pred)){

    if (is.element("curve",inputs)==TRUE) wCurve <- which(Curve$id==id.pred[i])
    if (is.element("scalar",inputs)==TRUE) wScalar <- which(Scalar$id==id.pred[i])
    if (is.element("factor",inputs)==TRUE) wFactor <- which(Factor$id==id.pred[i])

    noeud_courant <- 1

    while (is.element(noeud_courant, tree$feuilles)==FALSE){

      X <- get(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      type <- str_to_lower(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      var.split <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),3]))
      var.split.sum <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),4]))
      threshold <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),5]))

      # Maintenant il nous faut regarder la différence entre la moyenne à gauche et a droite et conclure :

      meanG <- tree$hist_nodes[[2*noeud_courant]]
      meanD <- tree$hist_nodes[[2*noeud_courant+1]]

      if (type=="curve"){

        fixed_var <- all.vars(X$model[[var.split]]$fixed)
        random_var <- all.vars(X$model[[var.split]]$random)
        model_var <- unique(c(fixed_var,random_var))

        data_model <- data.frame(id = as.numeric(X$id[wCurve]), time = X$time[wCurve],
                                 X$X[wCurve, , drop = FALSE])

        data_model <- data_model[,c("id",model_var)]

        RE <- predRE(tree$model_param[[noeud_courant]][[1]],
                     X$model[[var.split]], data_model)$bi

        ######################

        # autres resumes

        #####################

        data_summaries <- RE

        if (is.na(data_summaries[,var.split.sum])){
          noeud_courant <- NA
          break
        }

        if (data_summaries[,var.split.sum] < threshold){
          distG <- 0
          distD <- 1
        }else{
          distG <- 1
          distD <- 0
        }

      }
      if (type=="scalar"){

        if (X$X[wScalar,var.split] < threshold){
          distG <- 0
          distD <- 1
        }else{
          distG <- 1
          distD <- 0
        }

      }

      if (type=="factor"){
        distG <- -1*(is.element(X$X[wFactor,var.split],meanG))
        distD <- -1*(is.element(X$X[wFactor,var.split],meanD))
      }

      if (distG <= distD) { noeud_courant <- 2*noeud_courant}
      if (distD < distG) {noeud_courant <- 2*noeud_courant +1}


    }

    if(tree$Y$type=="curve" || tree$Y$type=="surv"){
      pred[i] <- noeud_courant
    }

    else{
      if(!is.na(noeud_courant)){
        pred[i] <- tree$Y_pred[[noeud_courant]]
      }else{
        pred[i] <- NA
      }
    }
  }
  return(pred)
}
