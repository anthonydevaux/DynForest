#' Predict with Dynamic random forests
#'
#' @param object : Dynamic random forest
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param aligned.shape [logical]:
#' @param timeScale [numeric]:
#' @param d_out [numeric]:
#' @param ... : optional parameters to be passed to the low level function
#'
#' @import kmlShape
#' @import stringr
#' @import RiemBase
#' @import Evomorph
#' @import geomorph
#'
#' @return
#' @export
#'
predict.DynForest <- function(object, Curve=NULL,Scalar=NULL,Factor=NULL,Shape=NULL, Image=NULL,aligned.shape=FALSE, timeScale=0.1, d_out=0.1,...){
  # La première étape est de toujours lire les prédicteurs ::

  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time,
                  model=object$curve.model)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }
  if (is.null(Shape)==FALSE){
    Shape <- list(type="shape",X=Shape$X,id=Shape$id)
  }
  if (is.null(Image)==FALSE){
    Image <- list(type="image",X=Image$X,id=Image$id)
  }

  ## Puis on prend les prédicteurs:

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs
  # On va les lires en mettant la maj sur les différents éléments qui le constituent :

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  if (is.element("shape",inputs)==TRUE & aligned.shape==FALSE){
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k],print.progress = FALSE)$coords
    }
    aligned.shape=TRUE
  }

  Id.pred <- unique(get(Inputs[1])$id)
  pred.feuille <- matrix(0, ncol(object$rf), length(Id.pred))

  if (object$type=="factor"){
    pred.feuille <- as.data.frame(matrix(0, ncol(object$rf), length(Id.pred)))
  }

  if (object$splitrule=="MM"){
    for (t in 1:ncol(object$rf)){
      pred.feuille[t,] <- pred.MMT(object$rf[,t], Curve = Curve,Scalar = Scalar,Factor=Factor,Shape=Shape,Image=Image, timeScale, aligned.shape = aligned.shape)
    }
  }

  if (object$type=="scalar"){
    pred <- apply(pred.feuille, 2, "mean")
    return(pred)
  }

  if (object$type=="factor"){
    pred.all <- apply(pred.feuille, 2, "table")
    val <- factor(rep(NA, length(pred.all)), levels=object$levels)
    proba <- rep(NA, length(pred.all))
    for (k in 1:length(pred.all)){
      val[k] <- factor(attributes(which.max(pred.all[[k]])))
      proba[k] <- max(pred.all[[k]])/ncol(object$rf)
    }
    prediction <- data.frame(pred=val, prob=proba)
    return(prediction)
  }

  if (object$type=="shape"){
    pred <- array(0, dim=c(object$size[1], object$size[2],length(Id.pred)))
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- array(0,dim=c(object$size[1],object$size[2],ncol(object$rf)))
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      Ms <- mshape(pred_courant)
      M <- matrix(0,dim(Ms)[1], dim(Ms)[2])
      M[,1] <- Ms[,1]
      M[,2] <- Ms[,2]
      pred[,,l] <- M
    }
  }

  if (object$type=="image"){
    pred <- array(0, dim=c(object$size[1], object$size[2],length(Id.pred)))
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- array(0,dim=c(object$size[1],object$size[2],ncol(object$rf)))
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant[,,k] <- object$rf[,k]$Y_pred[[pred.feuille[k,l]]]
      }
      donnees <- riemfactory(pred_courant[,,,drop=FALSE])
      pred[,,l] <- rbase.mean(donnees)$x
    }
  }

  if (object$type=="curve"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[,k]$Y_pred[[pred.feuille[k,l]]])[1]),object$rf[,k]$Y_pred[[pred.feuille[k,l]]]))
      }
      predouille <- kmlShape::meanFrechet(pred_courant, timeScale = timeScale)
      predouille <- cbind(predouille, rep(Id.pred[l],dim(predouille)[1]))
      pred <- rbind(pred, predouille)
    }
    names(pred) <- c("times", "traj", "ID")
  }

  if (object$type=="surv"){
    pred <- NULL
    for (l in 1:dim(pred.feuille)[2]){
      pred_courant <- NULL
      for(k in 1:dim(pred.feuille)[1]){
        pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[,k]$Y_pred[[pred.feuille[k,l]]])[1]),object$rf[,k]$Y_pred[[pred.feuille[k,l]]]))
      }
      Pred = cbind(sort(unique(pred_courant$times)), rep(NA,length(sort(unique(pred_courant$times)))))
      ### Maintenant il faut prédire à partir de cet élément::
      pred_courant2 = NULL
      id_courant = NULL
      for (k in sort(unique(pred_courant$times))){
        w= which(pred_courant$times >= k)
        for (j in unique(pred_courant[,1])){
          w_y = intersect(w,which(pred_courant[,1]==j))
          if (length(w_y)>= 1){
            pred_courant2 = c(pred_courant2, pred_courant[w_y,][which.min(pred_courant[w_y,2]),3])
          }
        }
        Pred[which(Pred[,1]==k),2] = mean(pred_courant2)
      }
      predouille <- cbind(Pred, rep(Id.pred[l],dim(Pred)[1]))
      pred <- rbind(pred, predouille)
    }
    #names(pred) <- c("times", "traj", "ID")
  }
  return(pred)
}
