#' Predict Mixed Model Tree
#'
#' @param tree : Mixed Model Tree
#' @param Curve [list]: A list that contains the input curves.
#' @param Scalar [list]: A list that contains the input scalars.
#' @param Factor [list]: A list that contains the input factors.
#' @param Shape [list]: A list that contains the input shape.
#' @param Image [list]: A list that contains the input images.
#' @param aligned.shape [logical]: \code{TRUE} if the input shapes are aligned and normalized (\code{aligned.shape=FALSE} by default)
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
pred.MMT <- function(tree, Curve=NULL,Scalar=NULL,Factor=NULL,Shape=NULL,Image=NULL, aligned.shape=FALSE ,timeScale=0.1){

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  id.pred <- unique(get(Inputs[1])$id)

  if (is.element("shape",inputs)==TRUE & aligned.shape==FALSE){
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k],print.progress = FALSE)$coords
    }
  }


  if (tree$Y$type=="factor"){
    pred <- factor(rep(NA, length(id.pred)),levels=tree$Ylevels)
  }

  else{pred <- rep(NA,length(id.pred))}

  for (i in 1:length(id.pred)){

    if (is.element("curve",inputs)==TRUE) wCurve <- which(Curve$id==id.pred[i])
    if (is.element("scalar",inputs)==TRUE) wScalar <- which(Scalar$id==id.pred[i])
    if (is.element("factor",inputs)==TRUE) wFactor <- which(Factor$id==id.pred[i])
    if (is.element("shape",inputs)==TRUE) wShape <- which(Shape$id==id.pred[i])
    if (is.element("image",inputs)==TRUE) wImage <- which(Image$id==id.pred[i])

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

        data_model <- data.frame(id = as.numeric(X$id[wCurve]), time = X$time[wCurve],
                                 X$X[wCurve, , drop = FALSE])

        RE <- predRE(tree$model_param[[noeud_courant]][[1]],
                     X$model[[var.split]], data_model)$bi

        ######################

        # autres resumes

        #####################

        data_summaries <- RE

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

      if (type=="shape"){
        elementz <- array(X$X[,,wShape,var.split],dim = c(nrow(meanG),ncol(meanG),1))
        distG <- ShapeDist(elementz,meanG)
        distD <- ShapeDist(elementz, meanD)
      }

      if (type=="image"){
        distG <- rbase.pdist2(riemfactory(array(data = meanG$x,dim=c(nrow(meanG$x), ncol(meanG$x),1))),riemfactory(array(X$X[,,wImage,var.split],dim=c(nrow(meanG$x), ncol(meanG$x),1))))
        distD <- rbase.pdist2(riemfactory(array(data = meanD$x,dim=c(nrow(meanD$x), ncol(meanD$x),1))),riemfactory(array(X$X[,,wImage,var.split],dim=c(nrow(meanD$x), ncol(meanD$x),1))))
      }

      if (type=="factor"){
        distG <- -1*(is.element(X$X[wFactor,var.split],meanG))
        distD <- -1*(is.element(X$X[wFactor,var.split],meanD))
      }

      if (distG <= distD) { noeud_courant <- 2*noeud_courant}
      if (distD < distG) {noeud_courant <- 2*noeud_courant +1}


    }

    if(tree$Y$type=="curve" || tree$Y$type=="image" || tree$Y$type=="shape" || tree$Y$type=="surv"){
      pred[i] <- noeud_courant
    }

    else{
      pred[i] <- tree$Y_pred[[noeud_courant]]
    }
  }
  return(pred)
}
