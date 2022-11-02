#' Plot results about the most predictive variables used in DynForest
#'
#' This function displays a plot of the most predictive variables with the minimal depth (for class \code{DynForestVarDepth}), the variable importance (for class \code{DynForestVIMP}) or the grouped variable importance (for class \code{DynForestgVIMP}).
#'
#' @param x Object inheriting from classes \code{DynForestVarDepth}, \code{DynForestVIMP} or \code{DynForestgVIMP}, to respectively plot the minimal depth, the variable importance or grouped variable importance.
#' @param plot_level For \code{DynForestVarDepth} object, compute the statistic at predictor (\code{plot_level}="predictor") or feature (\code{plot_level}="feature") level
#' @param PCT For \code{DynForestVIMP} or \code{DynForestgVIMP} object, display VIMP statistic in percentage. Default value is FALSE.
#' @param ordering For \code{DynForestVIMP} object, order predictors according to VIMP value. Default value is TRUE.
#' @param ... Optional parameters to be passed to the low level function
#'
#' @import ggplot2
#' @importFrom stringr str_order
#'
#' @seealso \code{\link{DynForest} \link{var_depth} \link{compute_VIMP} \link{compute_gVIMP}}
#'
#' @return \code{plot()} function displays: \tabular{ll}{
#'    With \code{DynForestVarDepth} \tab the minimal depth for each predictor/feature \cr
#'    \tab \cr
#'    With \code{DynForestVIMP} \tab the VIMP for each predictor \cr
#'    \tab \cr
#'    With \code{DynForestVarDepth} \tab the grouped-VIMP for each given group \cr
#' }
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
#'                      cause = 2, ncores = 1, seed = 1234)
#'
#' # Run var_depth function
#' res_varDepth <- var_depth(res_dyn)
#'
#' # Plot minimal depth
#' plot(x = res_varDepth, plot_level = "feature")
#'
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_VIMP(DynForest_obj = res_dyn, ncores = 2)
#'
#' # Plot VIMP
#' plot(x = res_dyn_VIMP, PCT = TRUE)
#'
#' # Compute gVIMP statistic
#' res_dyn_gVIMP <- compute_gVIMP(DynForest_obj = res_dyn,
#'                                group = list(group1 = c("serBilir","SGOT"),
#'                                             group2 = c("albumin","alkaline")),
#'                                ncores = 2)
#'
#' # Plot gVIMP
#' plot(x = res_dyn_gVIMP, PCT = TRUE)
#'
#' }
#'
#' @name plot.DynForest
#' @export
plot.DynForestVarDepth <- function(x, plot_level = c("predictor","feature"), ...){

  # checking
  if (!all(plot_level%in%c("predictor","feature"))){
    stop("Only 'predictor' and 'feature' options are allowed for plot_level argument!")
  }

  if (length(plot_level)>1){
    plot_level <- plot_level[1]
  }

  min_depth_all <- x$var_node_depth
  depth.df <- data.frame(var = rep(min_depth_all$var, ncol(min_depth_all)-1),
                         id_node = unlist(c(x$var_node_depth[,2:ncol(min_depth_all)])),
                         tree = rep(seq(ncol(min_depth_all)-1), each = nrow(min_depth_all)))
  depth.df$group <- sub("\\..*", "", depth.df$var)
  depth.df <- depth.df[order(depth.df$var),]

  if (plot_level=="feature"){

    depth.nbtree <- aggregate(id_node ~ var, data = depth.df, FUN = function(x){
      return(sum(!is.na(x)))
    })

    g <- ggplot(depth.df, aes_string(x = "var", y = "id_node")) +
      geom_boxplot(aes_string(fill = "group")) +
      geom_text(data = depth.nbtree, aes_string(x = "var", label = "id_node"),
                y = max(depth.df$id_node, na.rm = T) - 1) +
      scale_x_discrete(limits=rev(unique(depth.df$var)[str_order(unique(depth.df$var))])) +
      xlab("Features") +
      ylab("Minimal depth") +
      guides(fill = "none") +
      theme_bw() +
      theme(axis.title.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 14, face = "bold")) +
      coord_flip()

    return(print(g))

  }

  if (plot_level=="predictor"){

    depthVar.df <- aggregate(id_node ~ group + tree, data = depth.df, min, na.rm = T)

    depthVar.nbtree <- aggregate(id_node ~ group, data = depthVar.df, FUN = length)

    g <- ggplot(depthVar.df, aes_string(x = "group", y = "id_node")) +
      geom_boxplot(aes_string(fill = "group")) +
      geom_text(data = depthVar.nbtree, aes_string(x = "group", label = "id_node"),
                y = max(depthVar.df$id_node, na.rm = T) - 1) +
      scale_x_discrete(limits=rev(unique(depthVar.df$group)[str_order(unique(depthVar.df$group))])) +
      xlab("Predictors") +
      ylab("Minimal depth") +
      guides(fill = "none") +
      theme_bw() +
      theme(axis.title.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 14, face = "bold")) +
      coord_flip()

    return(print(g))

  }

}

#' @rdname plot.DynForest
#' @export
plot.DynForestVIMP <- function(x, PCT = FALSE, ordering = TRUE, ...){

  vimp.df <- data.frame(var = unlist(x$Inputs),
                        vimp = unlist(x$Importance))

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(x$tree_oob_err, na.rm = T) # vimp relative
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

#' @rdname plot.DynForest
#' @export
plot.DynForestgVIMP <- function(x, PCT = FALSE, ...){

  vimp.df <- data.frame(var = names(x$gVIMP),
                        vimp = x$gVIMP)

  if (PCT){
    vimp.df$vimp <- vimp.df$vimp*100/mean(x$tree_oob_err, na.rm = T) # vimp relative
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
