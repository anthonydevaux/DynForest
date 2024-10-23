#' Compute the importance of variables (VIMP) statistic
#'
#' @param dynforest_obj dynforest_obj \code{dynforest} object
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#' @param seed Seed to replicate results
#'
#' @importFrom methods is
#' @import doRNG
#'
#' @return \code{compute_vimp()} function returns a list with the following elements:\tabular{ll}{
#'    \code{Inputs} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Importance} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains a numeric vector of VIMP statistic predictor in \code{Inputs} value \cr
#'    \tab \cr
#'    \code{tree_oob_err} \tab A numeric vector containing the OOB error for each tree needed to compute the VIMP statistic \cr
#'    \tab \cr
#'    \code{IBS.range} \tab A vector containing the IBS min and max \cr
#' }
#'
#' @export
#'
#' @seealso [dynforest()]
#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Get Gaussian distribution for longitudinal predictors
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
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
#' # Run dynforest function
#' res_dyn <- dynforest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_vimp(dynforest_obj = res_dyn, ncores = 2, seed = 1234)
#'
#' }
compute_vimp <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                         ncores = NULL, seed = 1234){

  if (!methods::is(dynforest_obj,"dynforest")){
    cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You've supplied a {.cls {class(dynforest_obj)}} object"
    ))
  }

  if (dynforest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(dynforest_obj$data$Y$Y[,1])
    }
  }

  rf <- dynforest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  Y <- rf$data$Y
  timeVar <- rf$timeVar
  ntree <- ncol(rf$rf)
  Inputs <- names(rf$Inputs[!sapply(rf$Inputs,is.null)])

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
                     FUN=function(i){OOB.tree(rf$rf[,i], Longitudinal = Longitudinal, Numeric = Numeric, Factor = Factor, Y = Y,
                                              timeVar = timeVar, IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)},cl=cl)

  parallel::stopCluster(cl)

  # tree_oob_err <- rep(NA, ntree)
  # for (i in 1:ntree){
  #   tree_oob_err[i] = OOB.tree(rf$rf[,i], Longitudinal=Longitudinal,Numeric=Numeric,Factor = Factor, Y=Y,
  #                        IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)
  # }

  #####################

  Longitudinal.perm <- Longitudinal
  Numeric.perm <- Numeric
  Factor.perm <- Factor

  p <- NULL
  Importance.Longitudinal <- NULL
  Importance.Numeric <- NULL
  Importance.Factor <- NULL

  if (is.element("Longitudinal",Inputs)==TRUE){

    Longitudinal.err <- matrix(NA, ntree, ncol(Longitudinal$X))

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    Importance.Longitudinal <- foreach::foreach(p=1:ncol(Longitudinal$X),
                                         .combine = "c", .options.RNG = seed) %dorng% {
    # for (p in 1:ncol(Longitudinal$X)){

    Longitudinal.perm$X[,p] <- sample(x = na.omit(Longitudinal$X[,p]),
                                      size = length(Longitudinal$X[,p]),
                                      replace = TRUE) # avoid NA issue with permut

      for (k in 1:ntree){

        Longitudinal.err[k,p] <- OOB.tree(rf$rf[,k], Longitudinal = Longitudinal.perm, Numeric = Numeric, Factor = Factor, Y,
                                          timeVar = timeVar, IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

      }
      Longitudinal.perm$X[,p] <- Longitudinal$X[,p]
      res <- mean(Longitudinal.err[,p]- tree_oob_err)
    }

    parallel::stopCluster(cl)

  }


  if (is.element("Numeric",Inputs)==TRUE){

    Numeric.err <- matrix(NA, ntree, dim(Numeric$X)[2])

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    Importance.Numeric <- foreach::foreach(p=1:ncol(Numeric$X),
                                          .combine = "c", .options.RNG = seed) %dorng% {


    Numeric.perm$X[,p] <- sample(Numeric$X[,p])

    # for (p in 1:ncol(Numeric$X)){
      for (k in 1:ntree){

        Numeric.err[k,p] <- OOB.tree(rf$rf[,k], Longitudinal = Longitudinal, Numeric = Numeric.perm, Factor = Factor, Y,
                                     timeVar = timeVar, IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

      }
      Numeric.perm$X[,p] <- Numeric$X[,p]
      res <- mean(Numeric.err[,p]- tree_oob_err)
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
                                          .combine = "c", .options.RNG = seed) %dorng% {

    Factor.perm$X[,p] <- sample(Factor$X[,p])

    #for (p in 1:ncol(Factor$X)){

      for (k in 1:ntree){

        Factor.err[k,p] <- OOB.tree(rf$rf[,k], Longitudinal=Longitudinal, Numeric = Numeric, Factor=Factor.perm , Y,
                                    timeVar = timeVar, IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

      }

      Factor.perm$X[,p] <- Factor$X[,p]
      res <- mean(Factor.err[,p]- tree_oob_err)
    }

    parallel::stopCluster(cl)
  }

  Importance <- list(Longitudinal=as.vector(Importance.Longitudinal), Numeric=as.vector(Importance.Numeric), Factor=as.vector(Importance.Factor))

  out <- list(Inputs = dynforest_obj$Inputs,
              Importance = Importance,
              tree_oob_err = tree_oob_err,
              IBS.range = c(IBS.min, IBS.max))

  class(out) <- c("dynforestvimp")

  return(out)

}
