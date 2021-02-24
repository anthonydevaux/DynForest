#' OOB tree
#'
#' @param tree
#' @param Curve
#' @param Scalar
#' @param Factor
#' @param Shape
#' @param Image
#' @param Y
#' @param timeScale
#' @param d_out
#'
#' @import kmlShape
#' @import Evomorph
#' @import stringr
#' @import RiemBase
#'
#' @keywords internal
OOB.tree <- function(tree, Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL ,Y, timeScale=0.1, d_out=0.1,
                     splitrule=NULL){

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  BOOT <- tree$boot
  OOB <- setdiff(unique(Y$id), BOOT)
  xerror <- rep(NA,length(OOB))
  Scalar_courant <- NULL
  Factor_courant <- NULL
  Curve_courant <- NULL
  Image_courant <- NULL
  Shape_courant <- NULL

  if (Y$type=="curve" || Y$type=="surv"){
    for (i in OOB){
      id_wY <- which(Y$id== i)
      if (is.element("curve",inputs)==TRUE) {
        id_wXCurve <- which(Curve$id==i)
        Curve_courant <- list(type="curve",X=Curve$X[id_wXCurve,,drop=FALSE], id=Curve$id[id_wXCurve],time=Curve$time[id_wXCurve],
                              model=Curve$model)
      }

      if (is.element("shape",inputs)==TRUE){
        id_wXShape <- which(Shape$id==i)
        Shape_courant <- list(type="shape",X=Shape$X[,,id_wXShape,,drop=FALSE], id=Shape$id[id_wXShape])
      }

      if (is.element("image",inputs)==TRUE){
        id_wXImage <- which(Image$id==i)
        Image_courant <- list(type="image",X=Image$X[,,id_wXImage,,drop=FALSE], id=Image$id[id_wXImage])
      }

      if (is.element("factor",inputs)==TRUE){
        id_wXFactor <- which(Factor$id==i)
        Factor_courant <- list(type="factor",X=Factor$X[id_wXFactor,,drop=FALSE], id=Factor$id[id_wXFactor])
      }

      if (is.element("scalar",inputs)==TRUE){
        id_wXScalar <- which(Scalar$id==i)
        Scalar_courant <- list(type="scalar",X=Scalar$X[id_wXScalar,,drop=FALSE], id=Scalar$id[id_wXScalar])
      }

      if (splitrule=="MM"){
        pred_courant <- pred.MMT(tree, Curve=Curve_courant,Scalar=Scalar_courant,Factor=Factor_courant,Shape = Shape_courant,Image=Image_courant, timeScale=timeScale, aligned.shape = TRUE)
        #chancla <- DouglasPeuckerNbPoints(tree$Y_Curves[[pred_courant]]$times, tree$Y_Curves[[pred_courant]]$traj, nbPoints = length(stats::na.omit(Y[id_w])))
        xerror[which(OOB==i)] <- kmlShape::distFrechet(tree$Y_pred[[pred_courant]]$times, tree$Y_pred[[pred_courant]]$traj, Y$time[id_wY], Y$Y[id_wY], timeScale = d_out)^2
      }

    }
  }
  else {
    w_XCurve <- NULL
    w_XScalar <- NULL
    w_XFactor <- NULL
    w_XShape <- NULL
    w_XImage <- NULL
    w_y <- NULL
    for (i in OOB){

      if (is.element("curve",inputs)==TRUE) w_XCurve <- c(w_XCurve, which(Curve$id==i))
      if (is.element("scalar",inputs)==TRUE) w_XScalar <- c(w_XScalar, which(Scalar$id==i))
      if (is.element("factor",inputs)==TRUE) w_XFactor <- c(w_XFactor, which(Factor$id==i))
      if (is.element("shape",inputs)==TRUE) w_XShape <- c(w_XShape, which(Shape$id==i))
      if (is.element("image",inputs)==TRUE) w_XImage <- c(w_XImage, which(Image$id==i))

      w_y <- c(w_y, which(Y$id==i))
    }

    if (is.element("curve",inputs)==TRUE) Curve_courant <- list(type="curve",X=Curve$X[w_XCurve,,drop=FALSE], id=Curve$id[w_XCurve],time=Curve$time[w_XCurve],
                                                                model=Curve$model)
    if (is.element("scalar",inputs)==TRUE) Scalar_courant  <- list(type="scalar", X=Scalar$X[w_XScalar,, drop=FALSE], id=Scalar$id[w_XScalar])
    if (is.element("factor",inputs)==TRUE) Factor_courant  <- list(type="factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
    if (is.element("shape",inputs)==TRUE) Shape_courant  <- list(type="shape", X=Shape$X[,,w_XShape,, drop=FALSE], id=Shape$id[w_XShape])
    if (is.element("image",inputs)==TRUE) Image_courant  <- list(type="image", X=Image$X[,,w_XImage,, drop=FALSE], id=Image$id[w_XImage])

    pred <- pred.MMT(tree,Curve=Curve_courant,Scalar = Scalar_courant,Factor=Factor_courant, Shape=Shape_courant, Image = Image_courant, aligned.shape = TRUE)

    if (Y$type=="scalar"){xerror <- (Y$Y[w_y]-pred)^2}
    if (Y$type=="factor"){xerror <- 1*(pred!=Y$Y[w_y])}

    if (Y$type=="shape"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        xerror[l] <- ShapeDist(Y$Y[,,w_y[l], drop=FALSE],tree$Y_pred[[pred[l]]])^2
      }
    }

    if (Y$type=="image"){
      xerror <- rep(NA,length(pred))
      for (l in 1:length(pred)){
        moy_pred <- riemfactory(array(tree$Y_pred[[pred[l]]],dim=c(dim(Y$Y)[1],dim(Y$Y)[2],1)))
        vraie <- riemfactory(Y$Y[,,w_y[l], drop=FALSE])
        xerror[l] <- rbase.pdist2(vraie, moy_pred)^2
      }
    }

  }
  return(mean(xerror))
}
