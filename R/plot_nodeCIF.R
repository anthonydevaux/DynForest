#' Plot the estimated Cumulative Incidence Functions (CIF) for given tree nodes
#'
#' @param DynForest_obj A DynForest object from \code{DynForest()} function
#' @param tree Integer indicating the tree identifier
#' @param nodes Identifiers for the selected nodes
#'
#' @import ggplot2
#' @importFrom methods is
#'
#' @return Display the estimated CIF for given tree nodes
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
#' # Run DynForest function
#' res_dyn <- DynForest(timeData = timeData_train, fixedData = fixedData_train,
#'                      timeVar = "time", idVar = "id",
#'                      timeVarModel = timeVarModel, Y = Y,
#'                      ntree = 50, nodesize = 5, minsplit = 5,
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Display CIF for nodes 40 and 41 from tree 1
#' plot_nodeCIF(DynForest_obj = res_dyn, tree = 1, nodes = c(40,41))
#' }
#' @export
plot_nodeCIF <- function(DynForest_obj, tree = NULL, nodes = NULL){

  if (!methods::is(DynForest_obj,"DynForest")){
    stop("'DynForest_obj' should be a 'DynForestPred' class!")
  }

  if (!inherits(tree, "numeric")){
    stop("'tree' should be a numeric object containing the tree identifier!")
  }

  if (!all(inherits(nodes, "numeric"))){
    stop("'nodes' should be a numeric object containing the tree identifier!")
  }

  if (!any(tree==seq(DynForest_obj$param$ntree))){
    stop(paste0("'tree' should be chosen between 1 and ", DynForest_obj$param$ntree, "!"))
  }

  if (any(nodes>length(DynForest_obj$rf[,tree]$Y_pred))){
    stop("One selected node do not have CIF! Please verify the 'nodes' identifiers!")
  }

  if (any(sapply(nodes, FUN = function(x) is.null(DynForest_obj$rf[,tree]$Y_pred[[x]])))){
    stop("One selected node do not have CIF! Please verify the 'nodes' identifiers!")
  }

  # data transformation for ggplot2
  CIFs_nodes_list <- lapply(nodes, FUN = function(x){

    CIFs_node <- DynForest_obj$rf[,tree]$Y_pred[[x]]

    CIFs_node_list <- lapply(names(CIFs_node), FUN = function(y){

      CIF_node_cause <- CIFs_node[[y]]

      out <- data.frame(Node = rep(x, nrow(CIF_node_cause)),
                        Cause = rep(y, nrow(CIF_node_cause)),
                        Time = CIF_node_cause$times,
                        CIF = CIF_node_cause$traj)

      return(out)

    })

    return(do.call(rbind, CIFs_node_list))

  })

  data.CIF.plot <- do.call(rbind, CIFs_nodes_list)

  g <- ggplot(data.CIF.plot, aes_string(x = "Time", y = "CIF")) +
    geom_step(aes(group = Cause, color = Cause)) +
    facet_wrap(~ Node) +
    ylim(0,1) +
    theme_bw()

  print(g)

}
