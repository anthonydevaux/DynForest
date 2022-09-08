#' Compute the VIMP statistic
#'
#' @param DynForest_obj \code{DynForest} object
#' @param PCT Display VIMP statistic in pourcentage. Default value is FALSE.
#' @param ordering Order predictors according to VIMP value. Default value is TRUE.
#'
#' @import ggplot2
#'
#' @seealso \code{\link{DynForest} \link{compute_VIMP} \link{compute_gVIMP} \link{plot_gVIMP} \link{plot_mindepth}}
#'
#' @examples
#' \dontrun{
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
