#' Joint random forest for longitudinal and survival data
#'
#' DynForest function builds a joint random forest for survival analysis (Devaux et al., 2021), trajectory, regression or classification. Longitudinal data, factors and scalars are allowed as predictors.
#' Nodes are split on the candidate variable that maximize the splitting rule according to the outcome. To allow longitudinal data as candidate variable, mixed models are computed on those which are selected at each node.
#' Then, the random-effects are used are candidate variables. Out-of-bag prediction error and variable importance (VIMP) are also provided.
#'
#'
#' @param Curve [list]: A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar [list]: A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor [list]: A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Y [list]: A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{curve}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry [numeric]: Number of variables randomly drown as candidates at each split. The default value \code{p/3} but has to be tuned according to the OOB prediction error.
#' @param nodesize [numeric]
#' @param ntree [numeric]: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
#' @param ncores [numeric]: Number of cores used to build Frechet randomized trees in parallel, defaulting to number of cores of the computer minus 1.
#' @param timeScale [numeric]: Allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping. Only used when there are trajectories either in input or output.
#' @param imp [logical]: TRUE to compute the variables importance FALSE otherwise (default \code{imp=}TRUE)
#' @param d_out [string]: "euc" or "frec".
#' @param nsplit_option
#' @param cause
#' @param IBS.min
#' @param IBS.max
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
DynForest <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Y, mtry=NULL, ntree=100, ncores=NULL, timeScale=0.1, imp=TRUE, d_out=0.1,
                      nsplit_option = "quantile", nodesize = 1, cause = 1, IBS.min = 0, IBS.max = NULL, ...){


  if (Y$type=="surv"){
    Y$comp <- ifelse(length(unique(Y$Y[,2]))>2, TRUE, FALSE)
    causes <- sort(unique(Y$Y[which(Y$Y[,2]!=0),2]))
  }

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

  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  # On récupère le nombre de variables au total :
  nvar <- 0
  for (k in Inputs){
    nvar <- nvar + dim(get(k)$X)[length(dim(get(k)$X))]
  }

  if (is.null(mtry)==TRUE || mtry> nvar){
    mtry <- floor(nvar/3)*(floor(nvar/3)>=1) + 1*(floor(nvar/3)<1)
  }

  if(is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  print("Building Dynamic trees...")

  debut <- Sys.time()
  rf <-  rf_shape_para(Curve=Curve,Scalar=Scalar, Factor=Factor, Y=Y, mtry=mtry, ntree=ntree, timeScale = timeScale,ncores=ncores,
                       nsplit_option = nsplit_option, nodesize = nodesize, cause = cause)

  rf <- list(type=Y$type, rf=rf, levels=levels(Y$Y))

  print("OOB trees error...")
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  xerror <- pbsapply(1:ntree, FUN=function(i){OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y, timeScale=timeScale,
                                                       d_out=d_out, IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)},cl=cl)

  parallel::stopCluster(cl)

  # xerror <- rep(NA, ntree)
  # for (i in 1:ntree){
  #   cat(paste0("Tree ", i,"\n"))
  #   xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y, timeScale=timeScale, d_out=d_out,
  #                        IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)
  # }

  cat("OOB Forest error...")
  oob.err <- OOB.rfshape(rf,Curve = Curve,Scalar =Scalar,Factor=Factor, Y=Y, timeScale=timeScale, d_out=d_out,
                         IBS.min = IBS.min, IBS.max = IBS.max, cause = cause, ncores = ncores)
  cat("OK!\n")

  temps <- Sys.time() - debut

  cat("DynForest DONE!\n")

  # Ok pour le XERROR

  if (imp == FALSE && Y$type!="surv"){
    var.ini <- impurity(Y, timeScale)
    varex <- 1 - mean(oob.err$err)/var.ini
    drf <- list(rf=rf$rf,type=rf$type,levels=rf$levels, xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, varex=varex,
                Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
                Curve.model = Curve$model, comput.time=temps)
    class(drf) <- c("DynForest")
    return(drf)
  }

  if (imp == FALSE && Y$type=="surv"){
    drf <- list(rf=rf$rf,type=rf$type, times = sort(unique(c(0,Y$Y[,1]))), cause = cause, causes = causes,
                xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred,
                Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
                Curve.model = Curve$model, comput.time=temps)
    class(drf) <- c("DynForest")
    return(drf)
  }


  cat("Importance variables...\n")
  debut <- Sys.time()
  Curve.perm <- Curve
  Scalar.perm <- Scalar
  Factor.perm <- Factor

  Importance.Curve <- NULL
  Importance.Scalar <- NULL
  Importance.Factor <- NULL

  #X.perm <- list(type=X$type, X=X$X, id=X$id, time=X$time)
  if (is.element("curve",inputs)==TRUE){
    p=1
    cat("Curves...")
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


        Curve.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve.perm, Scalar = Scalar, Factor=Factor, Y, timeScale=timeScale,
                                   IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)

      }
      Curve.perm$X[,p] <- Curve$X[,p]
      res <- mean(Curve.err[,p]- xerror)
    }

    parallel::stopCluster(cl)

    cat("OK!\n")
  }


  if (is.element("scalar",inputs)==TRUE){
    p=1
    cat("Scalars...")
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

        Scalar.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar.perm, Factor=Factor, Y, timeScale=timeScale,
                                    IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)

      }
      Scalar.perm$X[,p] <- Scalar$X[,p]
      res <- mean(Scalar.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
    cat("OK!\n")
  }

  if (is.element("factor",inputs)==TRUE){
    p=1
    cat("Factors...")
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

        Factor.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor.perm , Y, timeScale=timeScale,
                                    IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)

      }
      ##on remet la variable en place :::
      Factor.perm$X[,p] <- Factor$X[,p]
      res <- mean(Factor.err[,p]- xerror)
    }

    parallel::stopCluster(cl)
    cat("OK!\n")
  }

  Importance <- list(Curve=as.vector(Importance.Curve), Scalar=as.vector(Importance.Scalar), Factor=as.vector(Importance.Factor))

  ############

  temps <- Sys.time() - debut

  cat("DynForest DONE!\n")

  if (Y$type == "surv"){
    drf <- list(rf=rf$rf,type=rf$type, times = sort(unique(c(0,Y$Y[,1]))), cause = cause, causes = causes,
                xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=Importance,
                Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
                Curve.model = Curve$model, comput.time=temps)
    class(drf) <- c("DynForest")
    return(drf)
  }
  var.ini <- impurity(Y, timeScale)
  varex <- 1 - mean(oob.err$err)/var.ini
  drf <- list(rf=rf$rf,type=rf$type,levels=rf$levels,xerror=xerror,oob.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=Importance, varex=varex,
              Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
              Curve.model = Curve$model, comput.time=temps)
  class(drf) <- c("DynForest")
  return(drf)
}
