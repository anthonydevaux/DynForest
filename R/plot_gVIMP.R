#' Title
#'
#' @param DynForest_obj
#' @param PCT
#'
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
plot_gVIMP <- function(DynForest_obj, PCT = FALSE){

  # checking
  if (!class(DynForest_obj)=="DynForest"){
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
