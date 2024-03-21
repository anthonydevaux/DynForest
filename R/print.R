#' Print function
#'
#' This function displays a brief summary regarding the trees (for class \code{DynForest}), a data frame with variable importance (for class \code{DynForestVIMP}) or the grouped variable importance (for class \code{DynForestgVIMP}).
#'
#' @param x Object inheriting from classes \code{DynForest}, \code{DynForestVIMP} or \code{DynForestgVIMP}.
#' @param ... Optional parameters to be passed to the low level function
#'
#'
#' @seealso \code{\link{DynForest} \link{var_depth} \link{compute_VIMP} \link{compute_gVIMP} \link{compute_OOBerror}}
#'
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
#' # Run DynForest function
#' res_dyn <- DynForest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Print function
#' print(res_dyn)
#'
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_VIMP(DynForest_obj = res_dyn, ncores = 2, seed = 1234)
#'
#' # Print function
#' print(res_dyn_VIMP)
#'
#' # Compute gVIMP statistic
#' res_dyn_gVIMP <- compute_gVIMP(DynForest_obj = res_dyn,
#'                                group = list(group1 = c("serBilir","SGOT"),
#'                                             group2 = c("albumin","alkaline")),
#'                                ncores = 2, seed = 1234)
#'
#' # Print function
#' print(res_dyn_gVIMP)
#'
#' # Run var_depth function
#' res_varDepth <- var_depth(res_dyn)
#'
#' # Print function
#' print(res_varDepth)
#'
#' }
#'
#' @rdname print.DynForest
#' @export
print.DynForest <- function(x, ...){

  if (!methods::is(x,"DynForest")){
    stop("'x' should be an object of 'DynForest' class!")
  }

  if (x$type=="surv"){
    if (length(x$causes)>1){
      type <- "survival (competing risk)"
    }else{
      type <- "survival"
      split.rule <- "Maximize log-rank statistic test"
    }
    oob.type <- "Integrated Brier Score"
    leaf.stat <- "Cumulative incidence function"
  }

  if (x$type=="factor"){
    type <- "classification"
  }

  if (x$type=="numeric"){
    type <- "continuous"
  }

  cat(paste0("DynForest in ", type, " mode"), "\n")
  cat("----------------","\n")
  cat(paste0("\t","Average depth per tree: ",
             round(mean(apply(x$rf, 2, FUN = function(x){
               mean(x$V_split$depth[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of leaves per tree: ",
             round(mean(apply(x$rf, 2, FUN = function(x){
               length(x$V_split$type[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of subjects per leaf: ",
             round(mean(apply(x$rf, 2, FUN = function(x){
               mean(x$V_split$N[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  if (type%in%c("survival (competing risk)","survival")){
    cat(paste0("\t","Average number of events of interest per leaf: ",
               round(mean(apply(x$rf, 2, FUN = function(x){
                 mean(x$V_split$Nevent[which(x$V_split$type=="Leaf")])}),
                 na.rm = T),2)),"\n")
  }
  cat("----------------","\n")
}

#' @rdname print.DynForest
#' @export
print.DynForestVIMP <- function(x, ...){

  if (!methods::is(x,"DynForestVIMP")){
    stop("'x' should be an object of 'DynForestVIMP' class!")
  }

  out <- data.frame("Predictors" = unlist(x$Inputs),
                    "Type" = rep(names(x$Inputs), sapply(x$Inputs, FUN = length)),
                    "VIMP" = unlist(x$Importance),
                    row.names = NULL)

  return(out)
}


#' @rdname print.DynForest
#' @export
print.DynForestgVIMP <- function(x, ...){

  if (!methods::is(x,"DynForestgVIMP")){
    stop("'x' should be an object of 'DynForestgVIMP' class!")
  }

  out <- data.frame("Group" = names(x$gVIMP),
                    "Predictors" = sapply(x$group, FUN = function(x) paste(x, collapse = " / ")),
                    "gVIMP" = x$gVIMP,
                    row.names = NULL)

  return(out)
}


#' @rdname print.DynForest
#' @export
print.DynForestVarDepth <- function(x, ...){

  if (!methods::is(x,"DynForestVarDepth")){
    stop("'x' should be an object of 'DynForestVarDepth' class!")
  }

  return(x$min_depth)
}


#' @rdname print.DynForest
#' @export
print.DynForestOOB <- function(x, ...){

  if (!methods::is(x,"DynForestOOB")){
    stop("'x' should be an object of 'DynForestOOB' class!")
  }

  return(mean(x$oob.err, na.rm = TRUE))
}


#' @rdname print.DynForest
#' @export
print.DynForestPred <- function(x, ...){

  if (!methods::is(x,"DynForestPred")){
    stop("'x' should be an object of 'DynForestPred' class!")
  }

  return(x$pred_indiv)
}
