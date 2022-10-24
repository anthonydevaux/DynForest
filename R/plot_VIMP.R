#' Compute the VIMP statistic
#'
#' @param DynForest_VIMP_obj \code{DynForest_VIMP} object
#' @param PCT Display VIMP statistic in percentage. Default value is FALSE.
#' @param ordering Order predictors according to VIMP value. Default value is TRUE.
#'
#' @import ggplot2
#'
#' @seealso \code{\link{DynForest} \link{compute_VIMP} \link{compute_gVIMP} \link{plot_gVIMP} \link{plot_mindepth}}
#'
#' @return Display the VIMP for each predictor
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
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_VIMP(DynForest_obj = res_dyn_OOB)
#'
#' # Plot VIMP
#' plot_VIMP(DynForest_VIMP_obj = res_dyn_VIMP)
#'
#' }
plot_VIMP <- function(DynForest_VIMP_obj, PCT = FALSE, ordering = TRUE){

  # checking
  if (!inherits(DynForest_VIMP_obj, "DynForest_VIMP")){
    stop("DynForest_VIMP_obj argument should be DynForest_VIMP class!")
  }

  vimp.df <- data.frame(var = unlist(DynForest_VIMP_obj$Inputs),
                        vimp = unlist(DynForest_VIMP_obj$Importance))

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(DynForest_VIMP_obj$tree_oob_err, na.rm = T) # vimp relative
  }

  if (ordering){
    g <- ggplot(vimp.df) +
      geom_bar(aes_string("var", "vimp"), stat = "identity") +
      scale_x_discrete(limits=vimp.df$var[order(vimp.df$vimp)]) +
      xlab("Predictors") +
      ylab(ifelse(PCT,"% VIMP","VIMP")) +
      coord_flip() +
      theme_bw()
  }else{
    g <- ggplot(vimp.df) +
      geom_bar(aes_string("var", "vimp"), stat = "identity") +
      xlab("Predictors") +
      ylab(ifelse(PCT,"% VIMP","VIMP")) +
      coord_flip() +
      theme_bw()
  }

  return(print(g))
}
