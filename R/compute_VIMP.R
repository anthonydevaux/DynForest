#' Compute the importance of variables (VIMP) statistic
#'
#' @param DynForest_obj \code{DynForest} object containing the dynamic random forest used on train data
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#'
#' @return \code{compute_OOBerror()} function return a list with the following elements:\tabular{ll}{
#'    \code{data} \tab A list containing the data used to grow the trees \cr
#'    \tab \cr
#'    \code{rf} \tab A table with each tree in column. Provide multiple charactistics about the tree building \cr
#'    \tab \cr
#'    \code{type} \tab Outcome type \cr
#'    \tab \cr
#'    \code{times} \tab A numeric vector containing the time-to-event for all subjects \cr
#'    \tab \cr
#'    \code{cause} \tab Indicating the cause of interest \cr
#'    \tab \cr
#'    \code{causes} \tab A numeric vector containing the causes indicator \cr
#'    \tab \cr
#'    \code{Inputs} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Curve.model} \tab A list of longitudinal markers containing the formula used for modeling in the random forest \cr
#'    \tab \cr
#'    \code{param} \tab A list containing the hyperparameters \cr
#'    \tab \cr
#'    \code{xerror} \tab A numeric vector containing the OOB error for each tree \cr
#'    \tab \cr
#'    \code{oob.err} \tab A numeric vector containing the OOB error for each subject \cr
#'    \tab \cr
#'    \code{oob.pred} \tab Outcome prediction for all subjects \cr
#'    \tab \cr
#'    \code{IBS.range} \tab A vector containing the IBS min and max \cr
#'    \tab \cr
#'    \code{Importance} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains a numeric vector of VIMP statistic predictor in \code{Inputs} value \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#' @export
#'
#' @examples
#' \dontrun{
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_VIMP(DynForest_obj = res_dyn_OOB)
#' }
compute_VIMP <- function(DynForest_obj, ncores = NULL){

  if (class(DynForest_obj)!="DynForest"){
    stop("'DynForest_obj' should be a 'DynForest' class!")
  }

  if (is.null(DynForest_obj$xerror)){
    stop("OOB error should be first computed using 'compute_OOBerror()' function!")
  }

  if (DynForest_obj$type=="surv"){
    IBS.min <- DynForest_obj$IBS.range[1]
    IBS.max <- DynForest_obj$IBS.range[2]
  }

  rf <- DynForest_obj
  Curve <- rf$data$Curve
  Scalar <- rf$data$Scalar
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  ntree <- ncol(rf$rf)
  Inputs <- names(rf$Inputs)

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  #####################

  cat("VIMP...")
  Curve.perm <- Curve
  Scalar.perm <- Scalar
  Factor.perm <- Factor

  Importance.Curve <- NULL
  Importance.Scalar <- NULL
  Importance.Factor <- NULL

  if (is.element("Curve",Inputs)==TRUE){

    cat("Curves...")
    Curve.err <- matrix(NA, ntree, ncol(Curve$X))

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    Importance.Curve <- foreach::foreach(p=1:ncol(Curve$X),
                                         #.packages = "kmlShape" ,
                                         .combine = "c") %dopar% {
      # for (p in 1:ncol(Curve$X)){

      for (k in 1:ntree){

        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)
        id_boot_Curve <- which(Curve$id%in%BOOT)

        # Il faut maintenant faire la permutation :

        Curve.perm$X[-id_boot_Curve,p] <- sample(x = na.omit(Curve.perm$X[-id_boot_Curve,p]),
                                                 size = length(Curve.perm$X[-id_boot_Curve,p]),
                                                 replace = TRUE) # avoid NA issue with permut

        Curve.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve.perm, Scalar = Scalar, Factor=Factor, Y,
                                   IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

      }
      Curve.perm$X[,p] <- Curve$X[,p]
      res <- mean(Curve.err[,p]- rf$xerror)
    }

    parallel::stopCluster(cl)

  }


  if (is.element("Scalar",Inputs)==TRUE){

    cat("Scalars...")
    Scalar.err <- matrix(NA, ntree, dim(Scalar$X)[2])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    Importance.Scalar <- foreach::foreach(p=1:ncol(Scalar$X),
                                          #.packages = "kmlShape" ,
                                          .combine = "c") %dopar% {

      for (k in 1:ntree){
        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)
        id_boot_Scalar <- which(Scalar$id%in%BOOT)

        Scalar.perm$X[-id_boot_Scalar,p] <- sample(Scalar.perm$X[-id_boot_Scalar,p])

        Scalar.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar.perm, Factor=Factor, Y,
                                    IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

      }
      Scalar.perm$X[,p] <- Scalar$X[,p]
      res <- mean(Scalar.err[,p]- rf$xerror)
    }

    parallel::stopCluster(cl)
  }

  if (is.element("Factor",Inputs)==TRUE){

    cat("Factors...")
    Factor.err <- matrix(NA, ntree, dim(Factor$X)[2])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    Importance.Factor <- foreach::foreach(p=1:ncol(Factor$X),
                                          #.packages = "kmlShape" ,
                                          .combine = "c") %dopar% {
    #for (p in 1:ncol(Factor$X)){

      for (k in 1:ntree){

        BOOT <- rf$rf[,k]$boot
        nboot <- length(unique(Y$id))- length(BOOT)
        id_boot_Factor <- which(Factor$id%in%BOOT)

        # Il faut maintenant faire la permutation :

        Factor.perm$X[-id_boot_Factor,p] <- sample(Factor.perm$X[-id_boot_Factor,p])

        Factor.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor.perm , Y,
                                    IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

      }
      ##on remet la variable en place :::
      Factor.perm$X[,p] <- Factor$X[,p]
      res <- mean(Factor.err[,p]- rf$xerror)
    }

    parallel::stopCluster(cl)
  }

  Importance <- list(Curve=as.vector(Importance.Curve), Scalar=as.vector(Importance.Scalar), Factor=as.vector(Importance.Factor))

  out <- list(rf = rf$rf, type = rf$type, times = rf$times, cause = rf$cause, causes = rf$causes,
              Inputs = rf$Inputs, Curve.model = rf$Curve.model, param = rf$param,
              comput.time = rf$comput.time,
              xerror = rf$xerror, oob.err = rf$oob.err, oob.pred = rf$oob.err,
              IBS.range = rf$IBS.range, Importance = Importance)

  cat("OK!\n")

  class(out) <- c("DynForest")

  return(out)

}
