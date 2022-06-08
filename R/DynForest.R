#' Random forest with multivariate longitudinal endogenous covariates
#'
#' Built a random forest using multivariate longitudinal endogenous covariates
#'
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{curve}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error. Default is `NULL`.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#' @param ntree Number of trees to grow. Default value set to 200.
#' @param ncores Number of cores used to grow trees in parallel. Default value is the number of cores of the computer-1.
#' @param OOB_error Compute the OOB error. Default value is TRUE.
#' @param imp Compute (1) the importance of variables (VIMP); (2) the importance of groups (gVIMP) if there is a \code{imp.group} argument. Default value is FALSE
#' @param imp.group A list of groups with the name of the predictors assigned in each group
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param IBS.min (Only with survival outcome) Minimal time to compute the Integrated Brier Score. Default value is set to 0.
#' @param IBS.max (Only with survival outcome) Maximal time to compute the Integrated Brier Score. Default value is set to the maximal time-to-event found.
#' @param seed Seed to replicate results
#' @param ... Optional parameters to be passed to the low level function
#'
#' @import stringr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import pbapply
#' @import stats
#' @import utils
#'
#' @details The function currently supports survival (competing or single event), continuous or factor outcome.
#'
#' FUTUR IMPROVEMENTS:
#'
#' @return DynForest function return a list with the following elements:\tabular{ll}{
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
#'    \code{xerror} \tab A numeric vector containing the OOB error for each tree \cr
#'    \tab \cr
#'    \code{oob.err} \tab A numeric vector containing the OOB error for each subject \cr
#'    \tab \cr
#'    \code{oob.pred} \tab Outcome prediction for all subjects \cr
#'    \tab \cr
#'    \code{Importance} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains a numeric vector of VIMP statistic predictor in \code{Inputs} value \cr
#'    \tab \cr
#'    \code{gVIMP} \tab A numeric vector containing the gVIMP for each group defined in \code{imp.group} argument \cr
#'    \tab \cr
#'    \code{Inputs} \tab A list of 3 elements: \code{Curve}, \code{Scalar} and \code{Factor}. Each element contains the names of the predictors \cr
#'    \tab \cr
#'    \code{Curve.model} \tab A list of longitudinal markers containing the formula used for modeling in the random forest \cr
#'    \tab \cr
#'    \code{comput.time} \tab Computation time \cr
#' }
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
#'
#' @references Devaux et al., Random survival forests for competing risks with multivariate longitudinal endogenous covariates (2022), arXiv
#'
#' @seealso \code{\link{predict.DynForest} \link{plot_VIMP} \link{plot_gVIMP} \link{plot_mindepth}}
#'
#' @examples
#' \dontrun{
#' data(pbc2)
#'
#' # Define time-independent continuous covariate
#' cont_covar <- list(X = pbc2_surv[,"age", drop = FALSE],
#'                   id = pbc2_surv$id)
#'
#' # Define time-independent non continuous covariates
#' fact_covar <- list(X = pbc2_surv[,c("drug","sex")],
#'                    id = pbc2_surv$id)
#'
#' # Define time-dependent continuous markers
#' cont_traj <- list(X = pbc2_long[,c("serBilir","serChol","albumin","alkaline")],
#'                   id = pbc2_long$id,
#'                   time = pbc2_long$time,
#'                   model = list(serBilir = list(fixed = serBilir ~ time,
#'                                                random = ~ time),
#'                                serChol = list(fixed = serChol ~ time + I(time^2),
#'                                               random = ~ time + I(time^2)),
#'                                albumin = list(fixed = albumin ~ time,
#'                                               random = ~ time),
#'                                alkaline = list(fixed = alkaline ~ time,
#'                                                random = ~ time))
#' )
#'
#' # Define outcome (survival here)
#' Y <- list(type = "surv",
#'           Y = Surv(pbc2_surv$years, factor(pbc2_surv$event)),
#'           id = pbc2_surv$id)
#'
#' # Run DynForest function
#' res_dyn <- DynForest(Curve = cont_traj, Factor = fact_covar, Scalar = cont_covar,
#'                      Y = Y, ntree = 200, imp = TRUE,
#'                      mtry = 4, nodesize = 2, minsplit = 3,
#'                      cause = 2)
#'
#' }
#' @export
#'
DynForest <- function(Curve=NULL,Scalar=NULL, Factor=NULL, Y, mtry=NULL, ntree=200, ncores=NULL,
                      OOB_error = TRUE, imp=FALSE, imp.group = NULL,
                      nsplit_option = "quantile", nodesize = 1, minsplit = 2, cause = 1,
                      IBS.min = 0, IBS.max = NULL, seed = round(runif(1,0,10000)),
                      ...){

  debut <- Sys.time()

  ###### Input checking ######
  if (is.null(Curve)==FALSE){
    Curve <- list(type="Curve",X=Curve$X,id=Curve$id,time=Curve$time,
                  model=Curve$model)
  }
  if (is.null(Scalar)==FALSE){
    Scalar <- list(type="Scalar",X=Scalar$X,id=Scalar$id)
  }
  if (is.null(Factor)==FALSE){
    Factor <- list(type="Factor",X=Factor$X,id=Factor$id)
  }

  Inputs <- read.Xarg(c(Curve,Scalar,Factor))
  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  ###### Output checking ######
  if (Y$type=="surv"){
    Y$comp <- ifelse(length(unique(Y$Y[,2]))>2, TRUE, FALSE)
    causes <- sort(unique(Y$Y[which(Y$Y[,2]!=0),2]))
  }

  ###### General checking ######
  if (!is.null(imp.group)){
    if (!all(unlist(imp.group)%in%c(colnames(Curve$X),colnames(Factor$X),colnames(Scalar$X)))){
      stop("Unknown predictor(s) from 'imp.group' argument!")
    }
  }

  Importance <- gVIMP <- NULL
  oob.err <- xerror <- NULL

  # number of predictors
  nvar <- sum(sapply(Inputs, FUN = function(x) ncol(get(x)$X)))

  # default mtry
  if (is.null(mtry)==TRUE || mtry> nvar){
    mtry <- floor(nvar/3)*(floor(nvar/3)>=1) + 1*(floor(nvar/3)<1)
  }

  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }

  ######### DynTree #########
  print("Growing Dynamic trees...")

  rf <-  rf_shape_para(Curve=Curve,Scalar=Scalar, Factor=Factor, Y=Y, mtry=mtry, ntree=ntree, ncores=ncores,
                       nsplit_option = nsplit_option, nodesize = nodesize, minsplit = minsplit, cause = cause, seed = seed)

  rf <- list(type=Y$type, rf=rf, levels=levels(Y$Y))

  ###########################

  if (OOB_error){

    print("OOB trees error...")
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    pck <- .packages()
    dir0 <- find.package()
    dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
    parallel::clusterExport(cl,list("pck","dir"),envir=environment())
    parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

    xerror <- pbsapply(1:ntree, FUN=function(i){OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y,
                                                         IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)},cl=cl)

    parallel::stopCluster(cl)

    # xerror <- rep(NA, ntree)
    # for (i in 1:ntree){
    #   cat(paste0("Tree ", i,"\n"))
    #   xerror[i] = OOB.tree(rf$rf[,i], Curve=Curve,Scalar=Scalar,Factor = Factor, Y=Y,
    #                        IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)
    # }

    cat("OOB Forest error...")
    oob.err <- OOB.rfshape(rf, Curve = Curve, Scalar = Scalar, Factor = Factor, Y = Y,
                           IBS.min = IBS.min, IBS.max = IBS.max, cause = cause, ncores = ncores)
    cat("OK!\n")

    if (imp){

      cat("Importance variables...\n")
      Curve.perm <- Curve
      Scalar.perm <- Scalar
      Factor.perm <- Factor

      Importance.Curve <- NULL
      Importance.Scalar <- NULL
      Importance.Factor <- NULL

      #X.perm <- list(type=X$type, X=X$X, id=X$id, time=X$time)
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
                                                 set.seed(123)
                                                 Curve.perm$X[-id_boot_Curve,p] <- sample(x = na.omit(Curve.perm$X[-id_boot_Curve,p]),
                                                                                          size = length(Curve.perm$X[-id_boot_Curve,p]),
                                                                                          replace = TRUE) # avoid NA issue with permut

                                                 Curve.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve.perm, Scalar = Scalar, Factor=Factor, Y,
                                                                            IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)

                                               }
                                               Curve.perm$X[,p] <- Curve$X[,p]
                                               res <- mean(Curve.err[,p]- xerror)
                                             }

        parallel::stopCluster(cl)

        cat("OK!\n")
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
                                                                              IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)

                                                }
                                                Scalar.perm$X[,p] <- Scalar$X[,p]
                                                res <- mean(Scalar.err[,p]- xerror)
                                              }

        parallel::stopCluster(cl)
        cat("OK!\n")
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

                                                for (k in 1:ntree){
                                                  BOOT <- rf$rf[,k]$boot
                                                  nboot <- length(unique(Y$id))- length(BOOT)
                                                  id_boot_Factor <- which(Factor$id%in%BOOT)

                                                  # Il faut maintenant faire la permutation :

                                                  Factor.perm$X[-id_boot_Factor,p] <- sample(Factor.perm$X[-id_boot_Factor,p])

                                                  Factor.err[k,p] <- OOB.tree(rf$rf[,k], Curve=Curve, Scalar = Scalar, Factor=Factor.perm , Y,
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
      # grouped VIMP

      if (!is.null(imp.group)){

        cat("Grouped VIMP...")

        gVIMP <- vector("numeric", length(imp.group))
        names(gVIMP) <- names(imp.group)

        for (g in 1:length(imp.group)){

          group <- imp.group[[g]]

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
                                                    IBS.min = IBS.min, IBS.max = IBS.max, cause = cause))

                                    # res[k] <- OOB.tree(rf$rf[,k], Curve = Curve.perm,
                                    #                    Scalar = Scalar.perm,
                                    #                    Factor = Factor.perm, Y,
                                    #                    IBS.min = IBS.min, IBS.max = IBS.max, cause = cause)

                                  }

          parallel::stopCluster(cl)

          gVIMP[g] <- mean(res - xerror)

        }

      }else{

        gVIMP <- NULL

      }

    }

  }

  cat("DynForest DONE!\n")

  if (Y$type == "surv"){
    drf <- list(rf = rf$rf, type = rf$type, times = sort(unique(c(0,Y$Y[,1]))), cause = cause, causes = causes,
                xerror = xerror, oob.err = oob.err$err, oob.pred = oob.err$oob.pred, Importance = Importance, gVIMP = gVIMP,
                Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
                Curve.model = Curve$model, param = list(mtry = mtry, nodesize = nodesize,
                                                        minsplit = minsplit, ntree = ntree),
                comput.time = Sys.time() - debut)
    class(drf) <- c("DynForest")
    return(drf)
  }
  var.ini <- impurity(Y)
  varex <- 1 - mean(oob.err$err)/var.ini
  drf <- list(rf = rf$rf, type = rf$type, levels = rf$levels, xerror = xerror, oob.err = oob.err$err, oob.pred = oob.err$oob.pred,
              Importance = Importance, gVIMP = gVIMP,
              varex = varex, Inputs = list(Curve = names(Curve$X), Scalar = names(Scalar$X), Factor = names(Factor$X)),
              Curve.model = Curve$model, param = list(mtry = mtry, nodesize = nodesize,
                                                      minsplit = NULL, ntree = ntree),
              comput.time = Sys.time() - debut)
  class(drf) <- c("DynForest")
  return(drf)
}
