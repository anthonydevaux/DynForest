#' Title
#'
#' @param DynForest_obj
#'
#' @importFrom stringr str_order
#' @return
#' @export
#'
#' @examples
var_depth <- function(DynForest_obj){

  if (any(DynForest_obj$Inputs%in%c("Curve"))){

    # mindepth and node at depth

    Curve_depth_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type=="Curve",])==0){
        return(data.frame("type" = character(), "num_noeud" = numeric()))
      }

      df <- aggregate(num_noeud ~ type + var_split + var_summary, x$V_split[x$V_split$type=="Curve",], min)
      df <- data.frame(type = paste(df$type, df$var_split, df$var_summary, sep = "."),
                       num_noeud = df$num_noeud)
      return(df)

    })

    Curve_depth <- aggregate(num_noeud ~ type,
                             do.call(rbind, Curve_depth_list), mean)
    Curve_depth <- Curve_depth[str_order(Curve_depth$type),]

    Curve_node_depth <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "type",
                                                          all.x = TRUE, all.y = TRUE),
                               Curve_depth_list)
    Curve_node_depth <- Curve_node_depth[str_order(Curve_node_depth$type),]
    colnames(Curve_node_depth)[2:ncol(Curve_node_depth)] <- paste0("tree", seq(ncol(Curve_node_depth)-1))

    # count by tree

    Curve_count_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type=="Curve",])==0){
        return(data.frame("type" = character(), "num_noeud" = numeric()))
      }

      df <- aggregate(num_noeud ~ type + var_split + var_summary, x$V_split[x$V_split$type=="Curve",], length)
      df <- data.frame(type = paste(df$type, df$var_split, df$var_summary, sep = "."),
                       count = df$num_noeud)
      return(df)

    })

    Curve_count <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "type",
                                                     all.x = TRUE, all.y = TRUE),
                          Curve_count_list)

    Curve_count <- Curve_count[str_order(Curve_count$type),]
    colnames(Curve_count)[2:ncol(Curve_count)] <- paste0("tree", seq(ncol(Curve_count)-1))

  }else{

    Curve_depth <- NULL
    Curve_node_depth <- NULL
    Curve_count <- NULL

  }

  if (any(DynForest_obj$Inputs%in%c("Scalar","Factor"))){

    # mindepth and node at depth

    Other_depth_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type!="Curve",])==0){
        return(data.frame("type" = character(), "num_noeud" = numeric()))
      }

      df <- aggregate(num_noeud ~ type + var_split, x$V_split[x$V_split$type!="Curve",], min)
      df <- data.frame(type = paste(df$type, df$var_split, sep = "."),
                       num_noeud = df$num_noeud)
      return(df)

    })

    Other_depth <- aggregate(num_noeud ~ type,
                             do.call(rbind, Other_depth_list), mean)
    Other_depth <- Other_depth[str_order(Other_depth$type),]

    Other_node_depth <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "type",
                                                          all.x = TRUE, all.y = TRUE),
                               Other_depth_list)
    Other_node_depth <- Other_node_depth[str_order(Other_node_depth$type),]
    colnames(Other_node_depth)[2:ncol(Other_node_depth)] <- paste0("tree", seq(ncol(Other_node_depth)-1))

    # count by tree

    Other_count_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type!="Curve",])==0){
        return(data.frame("type" = character(), "num_noeud" = numeric()))
      }

      df <- aggregate(num_noeud ~ type + var_split, x$V_split[x$V_split$type!="Curve",], length)
      df <- data.frame(type = paste(df$type, df$var_split, sep = "."),
                       count = df$num_noeud)
      return(df)

    })

    Other_count <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "type",
                                                     all.x = TRUE, all.y = TRUE),
                          Other_count_list)

    Other_count <- Other_count[str_order(Other_count$type),]
    colnames(Other_count)[2:ncol(Other_count)] <- paste0("tree", seq(ncol(Other_count)-1))

  }else{

    Other_depth <- NULL
    Other_node_depth <- NULL
    Other_count <- NULL

  }

  min_depth <- rbind(Curve_depth, Other_depth)
  min_depth$rank <- rank(min_depth$num_noeud)

  var_node_depth <- rbind(Curve_node_depth, Other_node_depth)
  rownames(var_node_depth) <- seq(nrow(var_node_depth))

  var_count <- rbind(Curve_count, Other_count)
  var_count[is.na(var_count)] <- 0
  rownames(var_count) <- seq(nrow(var_count))

  return(list(min_depth = min_depth, var_node_depth = var_node_depth,
              var_count = var_count))

}
