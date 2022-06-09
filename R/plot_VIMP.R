#' Compute the VIMP statistic
#'
#' @param DynForest_obj \code{DynForest} object
#' @param PCT Display VIMP statistic in pourcentage. Default value is FALSE.
#' @param ordering Order predictors according to VIMP value. Default value is TRUE.
#'
#' @import ggplot2
#'
#' @seealso \code{\link{DynForest} \link{plot_gVIMP} \link{plot_mindepth}}
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
#' # Plot VIMP
#' plot_VIMP(res_dyn)
#'
#' }
#' @export
plot_VIMP <- function(DynForest_obj, PCT = FALSE, ordering = TRUE){

  # checking
  if (!inherits(DynForest_obj, "DynForest")){
    stop("DynForest_obj argument should be DynForest class!")
  }

  vimp.df <- data.frame(var = unlist(DynForest_obj$Inputs),
                        vimp = unlist(DynForest_obj$Importance))

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(DynForest_obj$oob.err, na.rm = T) # vimp relative
  }

  if (ordering){
    g <- ggplot(vimp.df) +
      geom_bar(aes(reorder(var, vimp), vimp), stat = "identity") +
      xlab("Predictors") +
      ylab(ifelse(PCT,"% VIMP","VIMP")) +
      coord_flip() +
      theme_bw()
  }else{
    g <- ggplot(vimp.df) +
      geom_bar(aes(var, vimp), stat = "identity") +
      xlab("Predictors") +
      ylab(ifelse(PCT,"% VIMP","VIMP")) +
      coord_flip() +
      theme_bw()
  }

  return(print(g))
}
