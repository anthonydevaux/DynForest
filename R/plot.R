#' Plot function in dynforest
#'
#' This function displays a plot of CIF for a given node and tree (for class \code{dynforest}), the most predictive variables with the minimal depth (for class \code{dynforestvardepth}), the variable importance (for class \code{dynforestvimp}) or the grouped variable importance (for class \code{dynforestgvimp}).
#'
#' @param x Object inheriting from classes \code{dynforest}, \code{dynforestvardepth}, \code{dynforestvimp} or \code{dynforestgvimp}, to respectively plot the CIF, the minimal depth, the variable importance or grouped variable importance.
#' @param tree For \code{dynforest} class, integer indicating the tree identifier
#' @param nodes For \code{dynforest} class, identifiers for the selected nodes
#' @param id For \code{dynforest} and \code{dynforestpred} classes, identifier for a given subject
#' @param max_tree For \code{dynforest} class, integer indicating the number of tree to display while using \code{id} argument
#' @param plot_level For \code{dynforestvardepth} class, compute the statistic at predictor (\code{plot_level}="predictor") or feature (\code{plot_level}="feature") level
#' @param PCT For \code{dynforestvimp} or \code{dynforestgvimp} class, display VIMP statistic in percentage. Default value is FALSE.
#' @param ordering For \code{dynforestvimp} class, order predictors according to VIMP value. Default value is TRUE.
#' @param ... Optional parameters to be passed to the low level function
#'
#' @import ggplot2
#' @importFrom stringr str_order
#'
#' @seealso [dynforest()] [compute_ooberror()] [compute_vimp()] [compute_gvimp()] [compute_vardepth()]
#'
#' @return \code{plot()} function displays: \tabular{ll}{
#'    With \code{dynforestvardepth} \tab the minimal depth for each predictor/feature \cr
#'    \tab \cr
#'    With \code{dynforestvimp} \tab the VIMP for each predictor \cr
#'    \tab \cr
#'    With \code{dynforestgvimp} \tab the grouped-VIMP for each given group \cr
#' }
#'
#' @examples
#' \donttest{
#' data(pbc2)
#'
#' # Get Gaussian distribution for longitudinal predictors
#' pbc2$serBilir <- log(pbc2$serBilir)
#' pbc2$SGOT <- log(pbc2$SGOT)
#' pbc2$albumin <- log(pbc2$albumin)
#' pbc2$alkaline <- log(pbc2$alkaline)
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
#' # Run dynforest function
#' res_dyn <- dynforest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Plot estimated CIF at nodes 17 and 32
#' plot(x = res_dyn, tree = 1, nodes = c(17,32))
#'
#' # Run var_depth function
#' res_varDepth <- compute_vardepth(res_dyn)
#'
#' # Plot minimal depth
#' plot(x = res_varDepth, plot_level = "feature")
#'
#' # Compute VIMP statistic
#' res_dyn_VIMP <- compute_vimp(dynforest_obj = res_dyn, ncores = 2)
#'
#' # Plot VIMP
#' plot(x = res_dyn_VIMP, PCT = TRUE)
#'
#' # Compute gVIMP statistic
#' res_dyn_gVIMP <- compute_gvimp(dynforest_obj = res_dyn,
#'                                group = list(group1 = c("serBilir","SGOT"),
#'                                             group2 = c("albumin","alkaline")),
#'                                ncores = 2)
#'
#' # Plot gVIMP
#' plot(x = res_dyn_gVIMP, PCT = TRUE)
#'
#' # Sample 5 subjects to predict the event
#' set.seed(123)
#' id_pred <- sample(id, 5)
#'
#' # Create predictors objects
#' pbc2_pred <- pbc2[which(pbc2$id%in%id_pred),]
#' timeData_pred <- pbc2_pred[,c("id", "time", "serBilir", "SGOT", "albumin", "alkaline")]
#' fixedData_pred <- unique(pbc2_pred[,c("id","age","drug","sex")])
#'
#' # Predict the CIF function for the new subjects with landmark time at 4 years
#' pred_dyn <- predict(object = res_dyn,
#'                     timeData = timeData_pred, fixedData = fixedData_pred,
#'                     idVar = "id", timeVar = "time",
#'                     t0 = 4)
#'
#' # Plot predicted CIF for subjects 26 and 110
#' plot(x = pred_dyn, id = c(26, 110))
#'
#' }
#'
#' @rdname plot.dynforest
#' @export
plot.dynforest <- function(x, tree = NULL, nodes = NULL, id = NULL, max_tree = NULL, ...){

  if (!methods::is(x,"dynforest")){
    cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You've supplied a {.cls {class(dynforest_obj)}} object"
    ))
  }

  if (!is.null(tree)){

    if (!inherits(tree, "numeric")){
      cli_abort(c(
        "{.var tree} must be a numeric object containing the tree identifier",
        "x" = "You've supplied a {.cls {class(tree)}} object"
      ))
    }

    if (!any(tree==seq(x$param$ntree))){
      cli_abort(c(
        "{.var tree} must be chosen between 1 and {x$param$ntree}",
        "x" = "You've chosen {tree}"
      ))
    }

    if (all(!is.null(nodes))){
      if (!all(inherits(nodes, "numeric"))){
        cli_abort(c(
          "{.var nodes} must be a numeric vector containing the node identifiers",
          "x" = "You've supplied a {.cls {class(nodes)}} object"
        ))
      }
      if (!all(nodes%in%names(x$rf[,tree]$Y_pred))){
        cli_abort(c(
          "At least one selected node in {.var nodes} doesn't exist in {.var tree} {tree}"
        ))
      }
      if (any(sapply(nodes, FUN = function(node) is.null(x$rf[,tree]$Y_pred[[as.character(node)]])))){
        cli_abort(c(
          "At least one selected node in {.var nodes} doesn't exist in {.var tree} {tree}"
        ))
      }
    }else{
      nodes <- get_treenodes(dynforest_obj = x, tree = tree)
    }

    # data transformation for ggplot2
    CIFs_nodes_list <- lapply(nodes, FUN = function(node){

      CIFs_node <- x$rf[,tree]$Y_pred[[as.character(node)]]

      CIFs_node_list <- lapply(names(CIFs_node), FUN = function(y){

        CIF_node_cause <- CIFs_node[[y]]

        out <- data.frame(Node = rep(node, nrow(CIF_node_cause)),
                          Cause = rep(y, nrow(CIF_node_cause)),
                          Time = CIF_node_cause$times,
                          CIF = CIF_node_cause$traj)

        return(out)

      })

      return(do.call(rbind, CIFs_node_list))

    })

    data.CIF.plot <- do.call(rbind, CIFs_nodes_list)

    g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
      geom_step(aes_string(group = "Cause", color = "Cause")) +
      facet_wrap(~ Node) +
      ylim(0,1) +
      theme_bw()

    return(print(g))

  }

  if (!is.null(id)){

    nodes <- apply(x$rf, 2, FUN = function(y){
      leaf_tree <- y$leaves[which(y$idY==id)]
    })

    data.CIF.plot <- NULL

    for (tree_id in seq(length(nodes))){

      if (length(nodes[[tree_id]])>0){
        tree_node <- nodes[[tree_id]]
      }else{
        next()
      }

      CIFs_node <- x$rf[,tree_id]$Y_pred[[as.character(tree_node)]]

      CIFs_node_list <- lapply(names(CIFs_node), FUN = function(y){

        CIF_node_cause <- CIFs_node[[y]]

        out <- data.frame(Tree = rep(tree_id, nrow(CIF_node_cause)),
                          Node = rep(tree_node, nrow(CIF_node_cause)),
                          Cause = rep(y, nrow(CIF_node_cause)),
                          Time = CIF_node_cause$times,
                          CIF = CIF_node_cause$traj)

        return(out)

      })

      data.CIF.plot <- rbind(data.CIF.plot, do.call(rbind, CIFs_node_list))

    }

    data.CIF.plot$Tree_Node <- paste0("Tree ", data.CIF.plot$Tree, " / Node ", data.CIF.plot$Node)
    data.CIF.plot$Tree_Node <- factor(data.CIF.plot$Tree_Node, levels = unique(data.CIF.plot$Tree_Node))

    if (!is.null(max_tree)){
      max_tree_id <- unique(data.CIF.plot$Tree)[seq(max_tree)]
      data.CIF.plot <- data.CIF.plot[which(data.CIF.plot$Tree%in%max_tree_id),]
    }

    g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
      geom_step(aes_string(group = "Cause", color = "Cause")) +
      facet_wrap(~ Tree_Node) +
      ylim(0,1) +
      theme_bw()

    return(print(g))

  }

}


#' @name plot.dynforest
#' @export
plot.dynforestvardepth <- function(x, plot_level = c("predictor","feature"), ...){

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

#' @rdname plot.dynforest
#' @export
plot.dynforestvimp <- function(x, PCT = FALSE, ordering = TRUE, ...){

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

#' @rdname plot.dynforest
#' @export
plot.dynforestgvimp <- function(x, PCT = FALSE, ...){

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

#' @rdname plot.dynforest
#' @export
plot.dynforestpred <- function(x, id = NULL, ...){

  if (!methods::is(x,"dynforestpred")){
    stop("'x' should be a 'dynforestpred' class!")
  }

  if (is.null(id)){
    stop("'id' cannot be NULL!")
  }

  if (!all(id%in%rownames(x$pred_indiv))){
    stop("Predictions are not available for some subjects. Please verify the subjects identifiers!")
  }

  data.CIF <- x$pred_indiv
  data.CIF <- data.CIF[which(rownames(data.CIF)%in%id),, drop = FALSE]

  times <- x$times
  n.times <- length(times)

  data.CIF.plot <- data.frame(id = as.factor(rep(id, each = n.times)),
                              Time = rep(times, length(id)),
                              CIF = c(t(data.CIF)))

  g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
    geom_step(aes(group = id, color = id)) +
    ylim(0,1) +
    geom_vline(xintercept = x$t0, linetype = "dashed") +
    theme_bw()

  print(g)

}
