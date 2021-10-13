#' Function to update the list of parameters for each marker using those estimated from previous node
#'
#' @param current_node [integer] Current node of the tree
#' @param markers [string] Name of the markers to be updated
#' @param params [list] List to be updated for the current node
#'
#' @export
#' @return List with the updated parameters from the requested markers
#'
#' @examples
getParamMM <- function(current_node, markers, params){

  starting_node <- current_node
  params[[starting_node]] <- rep(list(NA), length(markers))
  names(params[[starting_node]]) <- markers

  while (current_node>1 & length(markers)>0){ # get initial values from upper nodes

    current_node <- current_node%/%2

    current_node_marker <-
      names(params[[current_node]])[which(names(params[[current_node]])%in%markers)]

    if (length(current_node_marker)>0){

      params[[starting_node]][current_node_marker] <- params[[current_node]][current_node_marker]
      markers <- markers[-which(markers%in%current_node_marker)]

    }

  }

  return(params)

}
