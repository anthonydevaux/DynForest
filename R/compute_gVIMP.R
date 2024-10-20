#' Compute the grouped importance of variables (gVIMP) statistic
#'
#' @inheritParams compute_vimp
#' @param group A list of groups with the name of the predictors assigned in each group
#'
#' @importFrom methods is
#' @import doRNG
#'
#' @return \code{compute_gvimp()} function returns a list with the following elements:\tabular{ll}{
#'    \code{Inputs} \tab A list of 3 elements: \code{Longitudinal}, \code{Numeric} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{group} \tab A list of each group defined in \code{group} argument \cr
#'    \tab \cr
#'    \code{gVIMP} \tab A numeric vector containing the gVIMP for each group defined in \code{group} argument \cr
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
#' # Compute gVIMP statistic
#' res_dyn_gVIMP <- compute_gvimp(dynforest_obj = res_dyn,
#'                                group = list(group1 = c("serBilir","SGOT"),
#'                                             group2 = c("albumin","alkaline")),
#'                                ncores = 2, seed = 1234)
#' }
compute_gvimp <- function(dynforest_obj, IBS.min = 0, IBS.max = NULL,
                          group = NULL, ncores = NULL, seed = 1234){

  if (!methods::is(dynforest_obj,"dynforest")){
    stop("'dynforest_obj' should be a 'dynforest' class!")
  }

  if (dynforest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(dynforest_obj$data$Y$Y[,1])
    }
  }

  if (is.null(group)){
    stop("'group' argument cannot be NULL! Please define groups to compute the gVIMP statistic!")
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
  #                        timeVar = timeVar, IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)
  # }

  #####################

  gVIMP <- vector("numeric", length(group))
  names(gVIMP) <- names(group)

  set.seed(seed) # set seed for permutation

  for (g in 1:length(group)){

    g_group <- group[[g]]

    Factor.perm <- Numeric.perm <- Longitudinal.perm <- NULL

    for (Input in Inputs){ # Duplicate Inputs for permutation

      assign(paste0(Input,".perm"), get(Input))

    }

    for (p in 1:length(g_group)){

      var_group <- g_group[p]

      # Factor permutation
      if (any(Inputs=="Factor")){
        if (any(var_group%in%colnames(Factor$X))){
          Factor.perm$X[, var_group] <- sample(Factor$X[, var_group])
        }
      }

      # Numeric permutation
      if (any(Inputs=="Numeric")){
        if (any(var_group%in%colnames(Numeric$X))){
          Numeric.perm$X[, var_group] <- sample(Numeric$X[, var_group])
        }
      }

      # Longitudinal permutation
      if (any(Inputs=="Longitudinal")){
        if (any(var_group%in%colnames(Longitudinal$X))){
          Longitudinal.perm$X[, var_group] <- sample(x = na.omit(Longitudinal$X[, var_group]),
                                                     size = length(Longitudinal$X[, var_group]),
                                                     replace = TRUE) # avoid NA issue with permut
        }
      }

    }

    k <- NULL

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    res <- foreach::foreach(k=1:ntree,
                            .combine = "c", .options.RNG = seed) %dorng% {

      # res <- vector("numeric", ntree)
      # for (k in 1:ntree){

      return(OOB.tree(rf$rf[,k], Longitudinal = Longitudinal.perm,
                      Numeric = Numeric.perm,
                      Factor = Factor.perm, Y,
                      timeVar = timeVar,
                      IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause))

      # res[k] <- OOB.tree(rf$rf[,k], Longitudinal = Longitudinal.perm,
      #                    Numeric = Numeric.perm,
      #                    Factor = Factor.perm, Y,
      #                    timeVar = timeVar,
      #                    IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

    }

    parallel::stopCluster(cl)

    gVIMP[g] <- mean(res - tree_oob_err)

  }

  out <- list(Inputs = dynforest_obj$Inputs,
              group = group,
              gVIMP = gVIMP,
              tree_oob_err = tree_oob_err,
              IBS.range = c(IBS.min, IBS.max))

  class(out) <- c("dynforestgvimp")

  return(out)
}
