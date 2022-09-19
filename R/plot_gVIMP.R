#' Compute the grouped-VIMP statistic
#'
#' @param DynForest_obj \code{DynForest} object
#' @param PCT Display VIMP statistic in pourcentage. Default value is FALSE.
#'
#' @import ggplot2
#'
#' @seealso \code{\link{DynForest} \link{compute_gVIMP} \link{compute_VIMP} \link{plot_VIMP} \link{plot_mindepth}}
#'
#' @return Display the grouped-VIMP for each given group
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
#'                                             group2 = c("albumin","alkaline")))
#'
#' # Plot gVIMP
#' plot_gVIMP(res_dyn_gVIMP)
#'
#' }
plot_gVIMP <- function(DynForest_obj, PCT = FALSE){

  # checking
  if (!inherits(DynForest_obj,"DynForest")){
    stop("DynForest_obj argument should be DynForest class!")
  }

  vimp.df <- data.frame(var = names(DynForest_obj$gVIMP),
                        vimp = DynForest_obj$gVIMP)

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(DynForest_obj$oob.err, na.rm = T) # vimp relative
  }

  g <- ggplot(vimp.df) +
    geom_bar(aes_string("var", "vimp"), stat = "identity") +
    scale_x_discrete(limits=vimp.df$var[order(vimp.df$vimp)]) +
    xlab("Group of predictors") +
    ylab(ifelse(PCT,"% grouped-VIMP","grouped-VIMP")) +
    coord_flip() +
    theme_bw()

  return(print(g))

}
