#' Compute the grouped importance of variables (gVIMP) statistic
#'
#' @param DynForest_obj \code{DynForest} object containing the dynamic random forest used on train data
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param group A list of groups with the name of the predictors assigned in each group
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#'
#' @importFrom methods is
#'
#' @return \code{compute_gVIMP()} function returns a list with the following elements:\tabular{ll}{
#'    \code{data} \tab A list containing the data used to grow the trees \cr
#'    \tab \cr
#'    \code{rf} \tab A table with each tree in column. Provide multiple characteristics about the tree building \cr
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
#'    \code{gVIMP} \tab A numeric vector containing the gVIMP for each group defined in \code{group} argument \cr
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
#' # Compute gVIMP statistic
#' res_dyn_gVIMP <- compute_gVIMP(DynForest_obj = res_dyn_OOB,
#'                                group = list(group1 = c("serBilir","SGOT"),
#'                                             group2 = c("albumin","alkaline")),
#'                                ncores = 2)
#' }
compute_gVIMP <- function(DynForest_obj, IBS.min = 0, IBS.max = NULL,
                          group = NULL, ncores = NULL){

  if (!methods::is(DynForest_obj,"DynForest")){
    stop("'DynForest_obj' should be a 'DynForest' class!")
  }

  if (DynForest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(DynForest_obj$data$Y$Y[,1])
    }
  }

  if (is.null(group)){
    stop("'group' argument cannot be NULL! Please define groups to compute the gVIMP statistic!")
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

  xerror <- pbsapply(1:ntree,
                     FUN=function(i){OOB.tree(rf$rf[,i], Curve = Curve, Scalar = Scalar, Factor = Factor, Y = Y,
                                              IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)},cl=cl)

  parallel::stopCluster(cl)

  # xerror <- rep(NA, ntree)
  # for (i in 1:ntree){
  #   xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y,
  #                        IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)
  # }

  #####################

  gVIMP <- vector("numeric", length(group))
  names(gVIMP) <- names(group)

  for (g in 1:length(group)){

    group <- group[[g]]

    id_boot_Curve <- id_boot_Factor <- id_boot_Scalar <- NULL
    Factor.perm <- Scalar.perm <- Curve.perm <- NULL

    for (Input in Inputs){ # id des id non bootstrap

      assign(paste0(Input,".perm"), get(Input))

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
                            .combine = "c") %dopar% {

      # res <- vector("numeric", ntree)
      # for (k in 1:ntree){

      BOOT <- rf$rf[,k]$boot
      nboot <- length(unique(Y$id))- length(BOOT)

      for (Input in Inputs){ # id des id non bootstrap

        assign(paste0("id_boot_", Input), which(get(Input)$id%in%BOOT))

      }

      for (p in 1:length(group)){

        var_group <- group[p]

        if (any(Inputs=="Factor")){
          if (any(var_group%in%colnames(Factor$X))){

            Factor.perm$X[-id_boot_Factor, var_group] <-
              sample(Factor.perm$X[-id_boot_Factor, var_group])

          }
        }

        if (any(Inputs=="Scalar")){
          if (any(var_group%in%colnames(Scalar$X))){

            Scalar.perm$X[-id_boot_Scalar, var_group] <-
              sample(Scalar.perm$X[-id_boot_Scalar, var_group])

          }
        }

        if (any(Inputs=="Curve")){
          if (any(var_group%in%colnames(Curve$X))){

            Curve.perm$X[-id_boot_Curve, var_group] <-
              sample(Curve.perm$X[-id_boot_Curve, var_group])

          }
        }

      }

      return(OOB.tree(rf$rf[,k], Curve = Curve.perm,
                      Scalar = Scalar.perm,
                      Factor = Factor.perm, Y,
                      IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause))

      # res[k] <- OOB.tree(rf$rf[,k], Curve = Curve.perm,
      #                    Scalar = Scalar.perm,
      #                    Factor = Factor.perm, Y,
      #                    IBS.min = IBS.min, IBS.max = IBS.max, cause = rf$cause)

    }

    parallel::stopCluster(cl)

    gVIMP[g] <- mean(res - rf$xerror)

  }

  out <- list(rf = rf$rf, type = rf$type, times = rf$times, cause = rf$cause, causes = rf$causes,
              Inputs = rf$Inputs, Curve.model = rf$Curve.model, param = rf$param,
              comput.time = rf$comput.time,
              oob.err = rf$oob.err, oob.pred = rf$oob.err,
              IBS.range = rf$IBS.range, gVIMP = gVIMP)

  class(out) <- c("DynForest")

  return(out)
}
