#' Compute the grouped importance of variables (gVIMP) statistic
#'
#' @param DynForest_obj \code{DynForest} object containing the dynamic random forest used on train data
#' @param group A list of groups with the name of the predictors assigned in each group
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
#'    \code{gVIMP} \tab A numeric vector containing the gVIMP for each group defined in \code{group} argument \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#' @export
#'
#' @examples
#' \dontrun{
#' # Compute gVIMP statistic
#' res_dyn_VIMP <- compute_gVIMP(DynForest_obj = res_dyn_OOB,
#'                               group = list(group1 = c("serBilir","SGOT"),
#'                                            group2 = c("albumin","alkaline")))
#' }
compute_gVIMP <- function(DynForest_obj, group = NULL, ncores = NULL){

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

  #####################

  cat("Grouped VIMP...")

  gVIMP <- vector("numeric", length(group))
  names(gVIMP) <- names(group)

  for (g in 1:length(group)){

    group <- group[[g]]

    id_boot_Curve <- id_boot_Factor <- id_boot_Scalar <- NULL
    Factor.perm <- Scalar.perm <- Curve.perm <- NULL

    for (Input in Inputs){ # id des id non bootstrap

      assign(paste0(Input,".perm"), get(Input))

    }

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
              xerror = rf$xerror, oob.err = rf$oob.err, oob.pred = rf$oob.err,
              IBS.range = rf$IBS.range, gVIMP = gVIMP)

  class(out) <- c("DynForest")

  cat("OK!\n")

  return(out)
}
