#' Display the summary of DynForest
#'
#' @param object \code{DynForest} or \code{DynForest_OOB} object
#' @param ... Optional parameters to be passed to the low level function
#'
#' @return Return some information about the random forest
#'
#' @author Anthony Devaux (\email{anthony.devaux@@u-bordeaux.fr})
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
#' # DynForest summary
#' summary(object = res_dyn_OOB)
#' }
#' @rdname summary.DynForest
#' @export
summary.DynForest <- function(object, ...){

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
    type <- "classification"
    oob.type <- "Missclassification"
    split.rule <- "Minimize variance"
    leaf.stat <- "Majority vote"
  }

  if (object$type=="scalar"){
    type <- "regression"
    oob.type <- "Mean square error"
    split.rule <- "Minimize Gini index"
    leaf.stat <- "Mean"
  }

  ##############################

  cat(paste0("DynForest executed with ", type, " mode"),"\n")
  cat(paste0("\t","Splitting rule: ", split.rule),"\n")
  cat(paste0("\t","Out-of-bag error type: ", oob.type),"\n")
  cat(paste0("\t","Leaf statistic: ", leaf.stat),"\n")
  cat("----------------","\n")

  # input
  cat("Input","\n")
  cat(paste0("\t","Number of subjects: ", length(unique(unlist(apply(object$rf, 2, FUN = function(x) x$idY))))),"\n")
  cat(paste0("\t","Curve: ", length(object$Inputs$Curve), " predictor(s)"),"\n")
  cat(paste0("\t","Scalar: ", length(object$Inputs$Scalar), " predictor(s)"),"\n")
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

  # DynForest summary
  cat("DynForest summary","\n")
  cat(paste0("\t","Average depth by tree: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               mean(x$V_split$depth[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of leaves by tree: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               length(x$V_split$type[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  cat(paste0("\t","Average number of subjects by leaf: ",
             round(mean(apply(object$rf, 2, FUN = function(x){
               mean(x$V_split$N[which(x$V_split$type=="Leaf")])}),
               na.rm = T),2)),"\n")
  if (type%in%c("survival (competing risk)","survival")){
    cat(paste0("\t","Average number of events of interest by leaf: ",
               round(mean(apply(object$rf, 2, FUN = function(x){
                 mean(x$V_split$Nevent[which(x$V_split$type=="Leaf")])}),
                 na.rm = T),2)),"\n")
  }

  cat("----------------","\n")

  # erreur out-of-bag
  cat(paste0("Out-of-bag error based on ", oob.type),"\n")
  cat(paste0("\t","Tree-based out-of-bag error: ",
             ifelse(!is.null(object$xerror),
                    round(mean(object$xerror, na.rm = T), 4),
                    "Not computed!")),"\n")
  cat(paste0("\t","Individual-based out-of-bag error: ",
             ifelse(!is.null(object$oob.err),
                    round(mean(object$oob.err, na.rm = T), 4),
                    "Not computed!")),"\n")
  cat("----------------","\n")

  # computation time
  cat("Time to build the random forest \n")
  cat("\t")
  print(object$comput.time)
  cat("\n")
  cat("----------------","\n")

}
