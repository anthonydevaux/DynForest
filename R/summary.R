#' Display the summary of dynforest
#'
#' @param object \code{dynforest} or \code{dynforestOOB} object
#' @param ... Optional parameters to be passed to the low level function
#'
#' @return Return some information about the random forest
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
#' # Compute OOB error
#' res_dyn_OOB <- compute_ooberror(dynforest_obj = res_dyn, ncores = 2)
#'
#' # dynforest summary
#' summary(object = res_dyn_OOB)
#' }
#' @rdname summary.dynforest
#' @export
summary.dynforest <- function(object, ...){

  if (object$type=="surv"){
    if (length(object$causes)>1){
      type <- "survival (competing risk)"
      split.rule <- "Fine & Gray statistic test"
    }else{
      type <- "survival"
      split.rule <- "Maximize log-rank statistic test"
    }
    oob.type <- "Integrated Brier Score"
    leaf.stat <- "Cumulative incidence function"
  }

  if (object$type=="factor"){
    type <- "categorical"
    oob.type <- "Missclassification"
    split.rule <- "Minimize weighted within-group Shannon entropy"
    leaf.stat <- "Majority vote"
  }

  if (object$type=="numeric"){
    type <- "continuous"
    oob.type <- "Mean square error"
    split.rule <- "Minimize weighted within-group variance"
    leaf.stat <- "Mean"
  }

  ##############################

  cat(paste0("dynforest executed for ", type, " outcome"),"\n")
  cat(paste0("\t","Splitting rule: ", split.rule),"\n")
  cat(paste0("\t","Out-of-bag error type: ", oob.type),"\n")
  cat(paste0("\t","Leaf statistic: ", leaf.stat),"\n")
  cat("----------------","\n")

  # input
  cat("Input","\n")
  cat(paste0("\t","Number of subjects: ", length(unique(unlist(apply(object$rf, 2, FUN = function(x) x$idY))))),"\n")
  cat(paste0("\t","Longitudinal: ", length(object$Inputs$Longitudinal), " predictor(s)"),"\n")
  cat(paste0("\t","Numeric: ", length(object$Inputs$Numeric), " predictor(s)"),"\n")
  cat(paste0("\t","Factor: ", length(object$Inputs$Factor), " predictor(s)"),"\n")
  cat("----------------","\n")

  # tuning parameters
  cat("Tuning parameters","\n")
  cat(paste0("\t","mtry: ", object$param$mtry),"\n")
  cat(paste0("\t","nodesize: ", object$param$nodesize),"\n")
  if (type%in%c("survival (competing risk)","survival")){
    cat(paste0("\t","minsplit: ", object$param$minsplit),"\n")
  }
  cat(paste0("\t","ntree: ", object$param$ntree),"\n")
  cat("----------------","\n")
  cat("----------------","\n")

  # dynforest summary
  cat("dynforest summary","\n")
  cat(paste0("\t","Average depth per tree: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               mean(x$V_split$depth[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of leaves per tree: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               length(x$V_split$type[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of subjects per leaf: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               mean(x$V_split$N[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  if (type%in%c("survival (competing risk)","survival")){
    cat(paste0("\t","Average number of events of interest per leaf: ",
               round(mean(apply(object$rf, 2, FUN = function(x){
                 mean(x$V_split$Nevent[which(x$V_split$type=="Leaf")])}),
                 na.rm = T),2)),"\n")
  }

  cat("----------------","\n")

  # computation time
  cat("Computation time \n")
  cat(paste0("\t","Number of cores used: ", object$ncores),"\n")
  cat("\t")
  print(object$comput.time)
  cat("----------------","\n")

}

#' @rdname summary.dynforest
#' @export
summary.dynforestoob <- function(object, ...){

  if (object$type=="surv"){
    if (length(object$causes)>1){
      type <- "survival (competing risk)"
      split.rule <- "Fine & Gray statistic test"
    }else{
      type <- "survival"
      split.rule <- "Maximize log-rank statistic test"
    }
    oob.type <- "Integrated Brier Score"
    leaf.stat <- "Cumulative incidence function"
  }

  if (object$type=="factor"){
    type <- "categorical"
    oob.type <- "Missclassification"
    split.rule <- "Minimize weighted within-group Shannon entropy"
    leaf.stat <- "Majority vote"
  }

  if (object$type=="numeric"){
    type <- "continuous"
    oob.type <- "Mean square error"
    split.rule <- "Minimize weighted within-group variance"
    leaf.stat <- "Mean"
  }

  ##############################

  cat(paste0("dynforest executed for ", type, " outcome"),"\n")
  cat(paste0("\t","Splitting rule: ", split.rule),"\n")
  cat(paste0("\t","Out-of-bag error type: ", oob.type),"\n")
  cat(paste0("\t","Leaf statistic: ", leaf.stat),"\n")
  cat("----------------","\n")

  # input
  cat("Input","\n")
  cat(paste0("\t","Number of subjects: ", length(unique(unlist(apply(object$rf, 2, FUN = function(x) x$idY))))),"\n")
  cat(paste0("\t","Longitudinal: ", length(object$Inputs$Longitudinal), " predictor(s)"),"\n")
  cat(paste0("\t","Numeric: ", length(object$Inputs$Numeric), " predictor(s)"),"\n")
  cat(paste0("\t","Factor: ", length(object$Inputs$Factor), " predictor(s)"),"\n")
  cat("----------------","\n")

  # tuning parameters
  cat("Tuning parameters","\n")
  cat(paste0("\t","mtry: ", object$param$mtry),"\n")
  cat(paste0("\t","nodesize: ", object$param$nodesize),"\n")
  if (type%in%c("survival (competing risk)","survival")){
    cat(paste0("\t","minsplit: ", object$param$minsplit),"\n")
  }
  cat(paste0("\t","ntree: ", object$param$ntree),"\n")
  cat("----------------","\n")
  cat("----------------","\n")

  # dynforest summary
  cat("dynforest summary","\n")
  cat(paste0("\t","Average depth per tree: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               mean(x$V_split$depth[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of leaves per tree: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               length(x$V_split$type[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of subjects per leaf: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               mean(x$V_split$N[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  if (type%in%c("survival (competing risk)","survival")){
    cat(paste0("\t","Average number of events of interest per leaf: ",
               round(mean(apply(object$rf, 2, FUN = function(x){
                 mean(x$V_split$Nevent[which(x$V_split$type=="Leaf")])}),
                 na.rm = T),2)),"\n")
  }

  cat("----------------","\n")

  # out-of-bag error
  cat(paste0("Out-of-bag error based on ", oob.type),"\n")
  cat(paste0("\t","Out-of-bag error: ", round(mean(object$oob.err, na.rm = T), 4)), "\n")
  cat("----------------","\n")

  # computation time
  cat("Computation time \n")
  cat(paste0("\t","Number of cores used: ", object$ncores),"\n")
  cat("\t")
  print(object$comput.time)
  cat("----------------","\n")

}
