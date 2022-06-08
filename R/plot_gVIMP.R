#' Compute the grouped-VIMP statistic
#'
#' @param DynForest_obj \code{DynForest} object
#' @param PCT Display VIMP statistic in pourcentage. Default value is FALSE.
#'
#' @import ggplot2
#'
#' @seealso \code{\link{DynForest} \link{plot_VIMP} \link{plot_mindepth}}
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
#'                      imp.group = list(group1 = c("serBilir", "serChol"),
#'                                       group2 = c("albumin", "alkaline")),
#'                      mtry = 4, nodesize = 2, minsplit = 3,
#'                      cause = 2)
#'
#' # Plot gVIMP
#' plot_gVIMP(res_dyn)
#'
#' }
#' @export
#'
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
    geom_bar(aes(reorder(var, vimp), vimp), stat = "identity") +
    xlab("Group of predictors") +
    ylab(ifelse(PCT,"% grouped-VIMP","grouped-VIMP")) +
    coord_flip() +
    theme_bw()

  return(print(g))

}
