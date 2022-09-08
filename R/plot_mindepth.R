#' Compute the average minimal depth statistic
#'
#' @param var_depth_obj Object from \code{var_depth} function
#' @param plot_level Compute the statistic at predictor (\code{plot_level}="predictor") or feature (\code{plot_level}="feature") level
#'
#' @import ggplot2
#' @importFrom stringr str_order
#'
#' @seealso \code{\link{DynForest} \link{var_depth} \link{plot_VIMP} \link{plot_gVIMP}}
#'
#' @examples
#' \dontrun{
#'
#' # Run var_depth function
#' res_varDepth <- var_depth(res_dyn)
#'
#' }
#'
#' @export
plot_mindepth <- function(var_depth_obj, plot_level = c("predictor","feature")){

  # checking
  if (!all(plot_level%in%c("predictor","feature"))){
    stop("Only 'predictor' and 'feature' options are allowed for plot_level argument!")
  }

  if (length(plot_level)>1){
    plot_level <- plot_level[1]
  }

  min_depth_all <- var_depth_obj$var_node_depth
  depth.df <- data.frame(var = rep(min_depth_all$var, ncol(min_depth_all)-1),
                         num_noeud = unlist(c(var_depth_obj$var_node_depth[,2:ncol(min_depth_all)])),
                         tree = rep(seq(ncol(min_depth_all)-1), each = nrow(min_depth_all)))
  depth.df$group <- sub("\\..*", "", depth.df$var)
  depth.df <- depth.df[order(depth.df$var),]

  if (plot_level=="feature"){

    depth.nbtree <- aggregate(num_noeud ~ var, data = depth.df, FUN = function(x){
      return(sum(!is.na(x)))
    })

    g <- ggplot(depth.df, aes_string(x = "var", y = "num_noeud")) +
      geom_boxplot(aes_string(fill = "group")) +
      geom_text(data = depth.nbtree, aes_string(x = "var", label = "num_noeud"),
                y = max(depth.df$num_noeud, na.rm = T) - 1) +
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

    depthVar.df <- aggregate(num_noeud ~ group + tree, data = depth.df, min, na.rm = T)

    depthVar.nbtree <- aggregate(num_noeud ~ group, data = depthVar.df, FUN = length)

    g <- ggplot(depthVar.df, aes_string(x = "group", y = "num_noeud")) +
      geom_boxplot(aes_string(fill = "group")) +
      geom_text(data = depthVar.nbtree, aes_string(x = "group", label = "num_noeud"),
                y = max(depthVar.df$num_noeud, na.rm = T) - 1) +
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
