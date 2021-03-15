#' Dynamic Random Forest
#'
#' This function builds Frechet random Forest introduced by Capitaine et.al, this includes the OOB predictions, OOB errors and variable importance computations.
#'
#'
#' @param Curve [list]: A list that contains the different input curves. It must contain the following elements (no choice): \code{X} the matrix of the different curves, each column code for a different curve variable; \code{id} is the vector of the identifiers for the different trajectories contained in \code{X}; \code{time} is the vector of the measurement times associated with the trajectories contained in \code{X}.
#' @param Scalar [list]: A list that contains the different input scalars. It must contain the following elements (no choice):  \code{X} the matrix of the scalars, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Factor [list]: A list that contains the different input factors. It must contain the following elements (no choice):  \code{X} the matrix of the factors, each column code for a different variable; \code{id} is the vector of the identifiers for each individual.
#' @param Shape [list]: A list that contains the different input shapes. It must contain the following elements (no choice):  \code{X} the array of the shapes of dimension \code{n}x2x\code{l}x\code{p} where \code{n} is the number of points for composing each shape, \code{l} is the number of shapes and \code{p} is the number of shapes variables, \code{id} is the vector of the identifiers for each individual.
#' @param Image [list]: A list that contains the different input images. It must contain the following elements (no choice):  \code{X} the array of the images of dimension \code{n}x\code{m}x\code{l}x\code{p} where \code{n}*\code{m} is the size of each image, \code{l} is the number of images and \code{p} is the number of shapes variables; \code{id} is the vector of the identifiers for each individual.
#' @param Y [list]: A list that contains the output, It must contain the following elements (no choice): \code{type} defines the nature of the output, can be "\code{curve}", "\code{sclalar}", "\code{factor}", "\code{shape}", "\code{image}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry [numeric]: Number of variables randomly sampled as candidates at each split. The default value \code{p/3}
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param ncores [numeric]: Number of cores used to build Frechet randomized trees in parallel, defaulting to number of cores of the computer minus 1.
#' @param timeScale [numeric]: Allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping. Only used when there are trajectories either in input or output.
#' @param imp [logical]: TRUE to compute the variables importance FALSE otherwise (default \code{imp=}TRUE)
#' @param d_out [string]: "euc" or "frec".
#' @param ... : optional parameters to be passed to the low level function
#'
#' @import stringr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import pbapply
#'
#' @return A Frechet random forest which is a list of the following elements: \itemize{
#' \item \code{rf:} a list of the \code{ntree} randomized Frechet trees that compose the forest.
#' \item \code{xerror :} a vector containing the OOB prediction error of each randomized Frechet tree composing the forest.
#' \item \code{OOB.err: } a vector containing the OOB prediction error of each individual in the learning sample.
#' \item \code{OOB.pred: } a list of the OOB prediction for each individual in the learning set.
#' \item \code{Importance: } A vector containing the variables importance.
#' \item \code{varex: } “pseudo R-squared”: Percentage of variance explained.
#' }
#' @export
#'
DynForest <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL ,Y, mtry=NULL, ntree=100,ncores=NULL, timeScale=0.1, imp=TRUE, d_out=0.1,
                      splitrule = "MM", nsplit_option = "quantile", nodesize = 1, ...){


  ### On va regarder les différentes entrées:
  if (is.null(Curve)==FALSE){
    Curve <- list(type="curve",X=Curve$X,id=Curve$id,time=Curve$time,
                  model=Curve$model)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="factor",X=Factor$X,id=Factor$id)
  }

  if (is.null(Shape)==FALSE){
    Shape <- list(type="shape",X=Shape$X,id=Shape$id)
    for (k in 1:dim(Shape$X)[length(dim(Shape$X))]){
      Shape$X[,,,k] <- gpagen(Shape$X[,,,k], print.progress = FALSE)$coords
    }
  }

  if (is.null(Image)==FALSE){
    Image <- list(type="image",X=Image$X,id=Image$id)
  }


  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }
  #

  if (Y$type=="shape"){
    Y$Y <- gpagen(Y$Y,print.progress = FALSE)$coords
  }

  # On récupère le nombre de variables au total :
  nvar <- 0
  for (k in Inputs){
    nvar <- nvar + dim(get(k)$X)[length(dim(get(k)$X))]
  }

  if (is.null(mtry)==TRUE || mtry> nvar){
    mtry <- floor(nvar/3)*(floor(nvar/3)>=1) + 1*(floor(nvar/3)<1)
  }

  if (is.null(Shape)!=TRUE || is.null(Image)!=TRUE) ERT <- TRUE

  size <- NULL
  if (Y$type=="shape" || Y$type=="image"){
    size <- c(dim(Y$Y)[1],dim(Y$Y)[2])
  }

  if(is.null(ncores)==TRUE){
    ncores <- detectCores()-1
  }

  print("Building the maximal Dynamic trees...")

  debut <- Sys.time()
  rf <-  rf_shape_para(Curve=Curve,Scalar=Scalar, Factor=Factor, Shape=Shape, Image=Image,Y=Y, mtry=mtry, ntree=ntree, timeScale = timeScale,ncores=ncores, aligned.shape = TRUE,
                       splitrule = splitrule, nsplit_option = nsplit_option, nodesize = nodesize)
  temps <- Sys.time() - debut

  if (Y$type=="shape" || Y$type=="image"){
    rf <- list(type=Y$type, rf=rf, size = dim(Y$Y) )
  }
  else {
    rf <- list(type=Y$type, rf=rf, levels=levels(Y$Y))
  }


  print("Forest constucted !")
  xerror <- rep(NA, ntree)

  print("calcul erreur oob")
  #cl <- parallel::makeCluster(ncores)
  #doParallel::registerDoParallel(cl)

  #xerror <- pbsapply(1:ntree, FUN=function(i){OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor,Shape=Shape,Image=Image, Y=Y_KM, timeScale=timeScale,d_out=d_out)},cl=cl)

  #parallel::stopCluster(cl)
  xerror = rep(NA,ntree)

  if (Y$type=="surv"){
    ykm = NULL
    tkm=NULL
    idkm = NULL

    for (i in unique(Y$id)){
      qui = which(Y$id==i)
      donnees = survfit(Surv(Y$time[qui],Y$Y[qui])~1)
      ykm=c(ykm, donnees$surv)
      tkm=c(tkm, donnees$time)
      idkm = c(idkm, rep(i,length(donnees$surv)))
    }
    Y= list(type="surv", Y=ykm, id=idkm, time=tkm)
  }

  for (i in 1:ntree){
    xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor,Shape=Shape,Image=Image, Y=Y, timeScale=timeScale,d_out=d_out,
                         splitrule = splitrule)
  }
  print("on passe erreur oob de la foret")
  oob.err <- OOB.rfshape(rf,Curve = Curve,Scalar =Scalar,Factor=Factor,Shape=Shape,Image=Image,Y=Y, timeScale=timeScale, d_out=d_out,
                         splitrule = splitrule)

  # Ok pour le XERROR

  if (imp == FALSE && Y$type!="surv"){
    var.ini <- impurity(Y, timeScale)
    varex <- 1 - mean(oob.err$err)/var.ini
    drf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, varex=varex, size=size, time=temps,
                splitrule=splitrule, curve.model=Curve$model, Inputs = Inputs)
    class(drf) <- c("DynForest")
    return(drf)
  }

  if (imp == FALSE && Y$type=="surv"){
    drf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, size=size, time=temps,
                splitrule=splitrule, curve.model=Curve$model, Inputs = Inputs)
    class(drf) <- c("DynForest")
    return(drf)
  }


  print("Importance calculation...")
  debut <- Sys.time()
  Curve.perm <- Curve
  Scalar.perm <- Scalar
  Factor.perm <- Factor
  Shape.perm <- Shape
  Image.perm <- Image

  Importance.Curve <- NULL
  Importance.Scalar <- NULL
  Importance.Factor <- NULL
  Importance.Shape <- NULL
  Importance.Image <- NULL

  #X.perm <- list(type=X$type, X=X$X, id=X$id, time=X$time)
  if (is.element("curve",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of curves')
    Curve.err <- matrix(NA, ntree, dim(Curve$X)[2])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    Importance.Curve <- foreach::foreach(p=1:dim(Curve$X)[2],.packages = "kmlShape" ,.combine = "c") %dopar% {
      for (k in 1:ntree){
        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Curve <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Curve <- c(id_boot_Curve, which(Curve$id==BOOT[i]))
        }

        # Il faut maintenant faire la permutation :

        Curve.perm$X[-id_boot_Curve,p] <- permutation_courbes(Curve$X[-id_boot_Curve,p], Curve$id[-id_boot_Curve])


        Curve.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve.perm, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale,
                                   splitrule = splitrule)

      }
      Curve.perm$X[,p] <- Curve$X[,p]
      res <- mean(Curve.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }


  if (is.element("scalar",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of scalars')
    Scalar.err <- matrix(NA, ntree, dim(Scalar$X)[2])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    Importance.Scalar <- foreach::foreach(p=1:dim(Scalar$X)[2],.packages = "kmlShape" ,.combine = "c") %dopar% {

      for (k in 1:ntree){
        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Scalar <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Scalar <- c(id_boot_Scalar, which(Scalar$id==BOOT[i]))
        }


        Scalar.perm$X[-id_boot_Scalar,p] <- sample(Scalar.perm$X[-id_boot_Scalar,p])

        Scalar.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar.perm, Factor=Factor,Shape=Shape, Image=Image, Y, timeScale=timeScale,
                                    splitrule = splitrule)

      }
      Scalar.perm$X[,p] <- Scalar$X[,p]
      res <- mean(Scalar.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }

  if (is.element("factor",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of factors')
    Factor.err <- matrix(NA, ntree, dim(Factor$X)[2])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    Importance.Factor <- foreach::foreach(p=1:dim(Factor$X)[2],.packages = "kmlShape" ,.combine = "c") %dopar% {

      for (k in 1:ntree){
        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Factor <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Factor <- c(id_boot_Factor, which(Factor$id==BOOT[i]))
        }

        # Il faut maintenant faire la permutation :

        Factor.perm$X[-id_boot_Factor,p] <- sample(Factor.perm$X[-id_boot_Factor,p])

        Factor.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor.perm ,Shape=Shape, Image=Image, Y, timeScale=timeScale,
                                    splitrule = splitrule)

      }
      ##on remet la variable en place :::
      Factor.perm$X[,p] <- Factor$X[,p]
      res <- mean(Factor.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }

  if (is.element("shape",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of shapes')
    Shape.err <- matrix(NA, ntree, dim(Shape$X)[length(dim(Shape$X))])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    Importance.Shape <- foreach::foreach(p=1:dim(Shape$X)[length(dim(Shape$X))],.packages = "kmlShape" ,.combine = "c") %dopar% {

      for (k in 1:ntree){
        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Shape <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Shape <- c(id_boot_Shape, which(Shape$id==BOOT[i]))
        }

        # Il faut maintenant faire la permutation :

        Shape.perm$X[,,-id_boot_Shape,p] <- permutation_shapes(Shape.perm$X[,,-id_boot_Shape, p], Shape.perm$id[-id_boot_Shape])

        Shape.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape.perm, Image=Image, Y, timeScale=timeScale,
                                   splitrule = splitrule)

      }
      ##on remet la variable en place :::
      Shape.perm$X[,,,p] <- Shape$X[,,,p]
      res <- mean(Shape.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }

  if (is.element("image",inputs)==TRUE){
    p=1
    print('Computing the importance on the space of images')
    Image.err <- matrix(NA, ntree, dim(Image$X)[length(dim(Image$X))])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    Importance.Image <- foreach::foreach(p=1:dim(Image$X)[length(dim(Image$X))],.packages = "kmlShape" ,.combine = "c") %dopar% {

      for (k in 1:ntree){
        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)

        id_boot_Image <- NULL
        for (i in 1:length(BOOT)){
          id_boot_Image <- c(id_boot_Image, which(Image$id==BOOT[i]))
        }

        # Il faut maintenant faire la permutation :

        Image.perm$X[,,-id_boot_Image,p] <- permutation_shapes(Image.perm$X[,,-id_boot_Image, p], Image.perm$id[-id_boot_Image])

        Image.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor,Shape=Shape, Image=Image.perm, Y, timeScale=timeScale,
                                   splitrule = splitrule)

      }
      ##on remet la variable en place :::
      Image.perm$X[,,,p] <- Image$X[,,,p]
      res <- mean(Image.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
  }

  Importance <- list(Curve=as.vector(Importance.Curve), Scalar=as.vector(Importance.Scalar), Factor=as.vector(Importance.Factor), Shape=as.vector(Importance.Shape), Image=as.vector(Importance.Image))

  ############

  temps.imp <- Sys.time() - debut

  if (Y$type == "surv"){
    drf <- list(rf=rf$rf,type=rf$type,levels=rf$levels,xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=Importance, time=temps, size=size,
                splitrule=splitrule, curve.model=Curve$model, Inputs = Inputs)
    class(drf) <- c("DynForest")
    return(drf)
  }
  var.ini <- impurity(Y, timeScale)
  varex <- 1 - mean(oob.err$err)/var.ini
  drf <- list(rf=rf$rf,type=rf$type,levels=rf$levels,xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=Importance, varex=varex, time=temps, size=size,
              splitrule=splitrule, curve.model=Curve$model, Inputs = Inputs)
  class(drf) <- c("DynForest")
  return(drf)
}
