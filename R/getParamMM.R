#' Function to update the list of parameters for each marker using those estimated from previous node
#'
#' @param current_node Current node of the tree
#' @param markers Character vector indicating the name of the markers to be updated
#' @param params List to be updated for the current node
#'
#' @return List with the updated parameters from the requested markers
#'
#' @keywords internal
#' @noRd
getParamMM <- function(current_node, markers, params){

  starting_node_chr <- as.character(current_node)
  params[[starting_node_chr]] <- rep(list(NA), length(markers))
  names(params[[starting_node_chr]]) <- markers

  while (current_node > 1 & length(markers) > 0) { # get initial values from upper nodes
    current_node <- current_node %/% 2
    current_node_chr <- as.character(current_node)

    current_node_marker_indices <- match(names(params[[current_node_chr]]), markers, nomatch = 0)
    current_node_marker <- names(params[[current_node_chr]])[current_node_marker_indices > 0]

    if (length(current_node_marker) > 0) {
      params[[starting_node_chr]][current_node_marker] <- params[[current_node_chr]][current_node_marker]
      markers <- setdiff(markers, current_node_marker)
    }
  }

  return(params)

}
