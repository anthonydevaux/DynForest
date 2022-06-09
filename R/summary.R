#' Display the summary of DynForest
#'
#' @param object \code{DynForest} object resulting from \code{DynForest} function
#' @param ... Optional parameters to be passed to the low level function
#' @return
#'
#' @examples
#' \dontrun{
#' data(pbc2)
#'
#' # Build survival data
#' pbc2_surv <- unique(pbc2[,c("id","age","drug","sex","years","event")])
#'
#' # Define time-independent continuous covariate
#' cont_covar <- list(X = pbc2_surv[,"age", drop = FALSE],
#'                    id = pbc2_surv$id)
#'
#' # Define time-independent non continuous covariates
#' fact_covar <- list(X = pbc2_surv[,c("drug","sex")],
#'                    id = pbc2_surv$id)
#'
#' # Define time-dependent continuous markers
#' cont_traj <- list(X = pbc2[,c("serBilir","SGOT","albumin","alkaline")],
#'                   id = pbc2$id,
#'                   time = pbc2$time,
#'                   model = list(serBilir = list(fixed = serBilir ~ time,
#'                                                random = ~ time),
#'                                SGOT = list(fixed = SGOT ~ time + I(time^2),
#'                                            random = ~ time + I(time^2)),
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
#'                      imp.group = list(group1 = c("serBilir","SGOT"),
#'                                       group2 = c("albumin","alkaline")),
#'                      mtry = 3, nodesize = 2, minsplit = 3,
#'                      cause = 2, seed = 1234)
#'
#' # DynForest summary
#'
#' summary(res_dyn)
#'
#' }
#'
#' @export
summary.DynForest <- function(object, ...){

  if (object$type=="surv"){
    if (length(object$causes)>1){
      type <- "survival (competing risk)"
      split.rule <- "Fine & Gray statistic test"
    }else{
      type <- "survival"
      split.rule <- "Log-rank statistic test"
    }
    oob.type <- "Integrated Brier Score"
    leaf.stat <- "Cumulative incidence function"
  }

  if (object$type=="factor"){
    type <- "classification"
    oob.type <- "XXXXXX"
    split.rule <- "XXXXXX"
  }

  if (object$type=="scalar"){
    type <- "regression"
    oob.type <- "XXXXXX"
    split.rule <- "XXXXX"
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
  cat(paste0("\t","Tree-based out-of-bag error: ", round(mean(object$xerror, na.rm = T), 4)),"\n")
  cat(paste0("\t","Individual-based out-of-bag error: ", round(mean(object$oob.err, na.rm = T), 4)),"\n")
  cat("----------------","\n")

  # computation time
  print(object$comput.time)
  cat("----------------","\n")

}
