#' Plot the individual Cumulative Incidence Function (CIF) for the interest cause
#'
#' @param DynForestPred_obj An DynForestPred object from \code{predict()} function
#' @param id Identifiers for the selected subjects
#'
#' @import ggplot2
#'
#' @return Display the CIF for selected subjects
#' @export
#'
#' @examples
#' \dontrun{
#' # Display CIF for subjects 102 and 260
#' plot_CIF(DynForestPred_obj = pred_dyn,
#'          id = c(102, 260))
#' }
plot_CIF <- function(DynForestPred_obj, id = NULL){

  if (class(DynForestPred_obj)!="DynForestPred"){
    stop("'DynForestPred_obj' should be a 'DynForestPred' class!")
  }

  if (is.null(id)){
    stop("'id' cannot be NULL!")
  }

  if (!all(id%in%rownames(DynForestPred_obj$pred_outcome))){
    stop("Predictions are not available for some subjects. Please verify the subjects identifiers!")
  }

  data.CIF <- DynForestPred_obj$pred_outcome
  data.CIF <- data.CIF[which(rownames(data.CIF)%in%id),, drop = FALSE]

  times <- DynForestPred_obj$times
  n.times <- length(times)

  data.CIF.plot <- data.frame(id = as.factor(rep(id, each = n.times)),
                              Time = rep(times, length(id)),
                              CIF = c(t(data.CIF)))

  g <- ggplot(data.CIF.plot, aes(x = Time, y = CIF)) +
    geom_step(aes(group = id, color = id)) +
    ylim(0,1) +
    geom_vline(xintercept = DynForestPred_obj$t0, linetype = "dashed") +
    theme_bw()

  print(g)

}
