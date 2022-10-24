#' Compute the importance of variables (VIMP) statistic
#'
#' @param DynForest_obj \code{DynForest} object containing the dynamic random forest used on train data
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#'
#' @importFrom methods is
#'
#' @return \code{compute_VIMP()} function returns a list with the following elements:\tabular{ll}{
#'    \code{Inputs} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Importance} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains a numeric vector of VIMP statistic predictor in \code{Inputs} value \cr
#'    \tab \cr
#'    \code{tree_oob_err} \tab A numeric vector containing the OOB error for each tree \cr
#'    \tab \cr
#'    \code{IBS.range} \tab A vector containing the IBS min and max \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#' @export
#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Sample 100 subjects
#' set.seed(1234)
#' id <- unique(pbc2$id)
#' id_sample <- sample(id, 100)
#' id_row <- which(pbc2$id%in%id_sample)
#'
#' pbc2_train <- pbc2[id_row,]
#'
#  Build longitudinal data
#' timeData_train <- pbc2_train[,c("id","time",
#'                                 "serBilir","SGOT",
#'                                 "albumin","alkaline")]
#'
#' # Create object with longitudinal association for each predictor
#' timeVarModel <- list(serBilir = list(fixed = serBilir ~ time,
#'                                      random = ~ time),
#'                      SGOT = list(fixed = SGOT ~ time + I(time^2),
#'                                  random = ~ time + I(time^2)),
#'                      albumin = list(fixed = albumin ~ time,
#'                                     random = ~ time),
#'                      alkaline = list(fixed = alkaline ~ time,
#'                                      random = ~ time))
#'
#' # Build fixed data
#' fixedData_train <- unique(pbc2_train[,c("id","age","drug","sex")])
#'
#' # Build outcome data
#' Y <- list(type = "surv",
#'           Y = unique(pbc2_train[,c("id","years","event")]))
#'
#' # Run DynForest function
#' res_dyn <- DynForest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Compute OOB error
#' res_dyn_OOB <- compute_OOBerror(DynForest_obj = res_dyn, ncores = 2)
#'
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_VIMP(DynForest_obj = res_dyn_OOB, ncores = 2)
#' }
compute_VIMP <- function(DynForest_obj, IBS.min = 0, IBS.max = NULL,
                         ncores = NULL){

  if (!methods::is(DynForest_obj,"DynForest")){
    stop("'DynForest_obj' should be a 'DynForest' class!")
  }

  if (DynForest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(DynForest_obj$data$Y$Y[,1])
    }
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

  ##############################

  pbapply::pboptions(type="none")

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  pck <- .packages()
  dir0 <- find.package()
  dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
  parallel::clusterExport(cl,list("pck","dir"),envir=environment())
  parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

  tree_oob_err <- pbsapply(1:ntree,
                     FUN=function(i){OOB.tree(rf$rf[,i], Curve = Curve, Scalar = Scalar, Factor = Factor, Y = Y,
                                              IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)},cl=cl)

  parallel::stopCluster(cl)

  # tree_oob_err <- rep(NA, ntree)
  # for (i in 1:ntree){
  #   tree_oob_err[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y,
  #                        IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)
  # }

  #####################

  Curve.perm <- Curve
  Scalar.perm <- Scalar
  Factor.perm <- Factor

  p <- NULL
  Importance.Curve <- NULL
  Importance.Scalar <- NULL
  Importance.Factor <- NULL

  if (is.element("Curve",Inputs)==TRUE){

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
      res <- mean(Curve.err[,p]- rf$tree_oob_err)
    }

    parallel::stopCluster(cl)

  }


  if (is.element("Scalar",Inputs)==TRUE){

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
      res <- mean(Scalar.err[,p]- rf$tree_oob_err)
    }

    parallel::stopCluster(cl)
  }

  if (is.element("Factor",Inputs)==TRUE){

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
      res <- mean(Factor.err[,p]- rf$tree_oob_err)
    }

    parallel::stopCluster(cl)
  }

  Importance <- list(Curve=as.vector(Importance.Curve), Scalar=as.vector(Importance.Scalar), Factor=as.vector(Importance.Factor))

  out <- list(Inputs = Inputs,
              Importance = Importance,
              tree_oob_err = tree_oob_err,
              IBS.range = c(IBS.min, IBS.max))

  class(out) <- c("DynForest_VIMP")

  return(out)

}
