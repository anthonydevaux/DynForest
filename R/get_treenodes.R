#' Extract nodes identifiers for a given tree
#'
#' @inheritParams get_tree
#'
#' @importFrom methods is
#'
#' @return Extract nodes identifiers for a given tree
#'
#' @seealso [dynforest()]
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
#' # Extract nodes identifiers for a given tree
#' get_treenodes(dynforest_obj = res_dyn, tree = 1)
#' }
#' @export
get_treenodes <- function(dynforest_obj, tree = NULL){

  if (!methods::is(dynforest_obj,"dynforest")){
    stop("'dynforest_obj' should be a 'dynforestPred' class!")
  }

  if (!inherits(tree, "numeric")){
    stop("'tree' should be a numeric object containing the tree identifier!")
  }

  if (!any(tree==seq(dynforest_obj$param$ntree))){
    stop(paste0("'tree' should be chosen between 1 and ", dynforest_obj$param$ntree, "!"))
  }

  tree_split <- get_tree(dynforest_obj = dynforest_obj, tree = tree)
  nodes_id <- tree_split$id_node[which(tree_split$type=="Leaf")]

  return(nodes_id)
}
