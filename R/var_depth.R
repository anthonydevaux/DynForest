#' Extract characteristics from the trees building process
#'
#' @param DynForest_obj \code{DynForest} object
#'
#' @importFrom stringr str_order
#'
#' @return var_depth function return a list with the following elements:\tabular{ll}{
#'    \code{min_depth} \tab A table providing for each feature in row: the average depth and the rank \cr
#'    \tab \cr
#'    \code{var_node_depth} \tab A table providing for each tree in column the minimal depth for each feature in row. NA indicates that the feature was not used for the corresponding tree \cr
#'    \tab \cr
#'    \code{var_count} \tab A table providing for each tree in column the number of times where the feature is used (in row). 0 value indicates that the feature was not used for the corresponding tree \cr
#' }
#'
#' @seealso \code{\link{DynForest}}
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
#'                      cause = 2, ncores = 2, seed = 1234)
#'
#' # Run var_depth function
#' res_varDepth <- var_depth(res_dyn)
#'
#' }
#'
#' @export
var_depth <- function(DynForest_obj){

  Inputs <- names(DynForest_obj$Inputs)[unlist(lapply(DynForest_obj$Inputs, FUN = function(x) return(!is.null(x))))]

  if (any(Inputs%in%c("Longitudinal"))){

    # mindepth and node at depth

    Longitudinal_depth_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type=="Longitudinal",])==0){
        return(data.frame("var" = character(), "depth" = numeric()))
      }

      df <- aggregate(depth ~ type + var_split + feature, x$V_split[x$V_split$type=="Longitudinal",], min)
      df <- data.frame(var = paste0(DynForest_obj$Inputs$Longitudinal[df$var_split], ".bi", df$feature-1),
                       depth = df$depth)
      return(df)

    })

    Longitudinal_depth <- aggregate(depth ~ var,
                             do.call(rbind, Longitudinal_depth_list), mean)
    Longitudinal_depth <- Longitudinal_depth[str_order(Longitudinal_depth$var),]

    Longitudinal_node_depth <- suppressWarnings(
      Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "var",
                                        all.x = TRUE, all.y = TRUE),
             Longitudinal_depth_list))
    Longitudinal_node_depth <- Longitudinal_node_depth[str_order(Longitudinal_node_depth$var),]
    colnames(Longitudinal_node_depth)[2:ncol(Longitudinal_node_depth)] <- paste0("tree", seq(ncol(Longitudinal_node_depth)-1))

    # count by tree

    Longitudinal_count_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type=="Longitudinal",])==0){
        return(data.frame("var" = character(), "depth" = numeric()))
      }

      df <- aggregate(depth ~ type + var_split + feature, x$V_split[x$V_split$type=="Longitudinal",], length)
      df <- data.frame(var = paste0(DynForest_obj$Inputs$Longitudinal[df$var_split], ".bi", df$feature-1),
                       depth = df$depth)
      return(df)

    })

    Longitudinal_count <- suppressWarnings(
      Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "var",
                                        all.x = TRUE, all.y = TRUE),
             Longitudinal_count_list))

    Longitudinal_count <- Longitudinal_count[str_order(Longitudinal_count$var),]
    colnames(Longitudinal_count)[2:ncol(Longitudinal_count)] <- paste0("tree", seq(ncol(Longitudinal_count)-1))

  }else{

    Longitudinal_depth <- NULL
    Longitudinal_node_depth <- NULL
    Longitudinal_count <- NULL

  }

  if (any(Inputs%in%c("Numeric","Factor"))){

    # mindepth and node at depth

    Other_depth_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type%in%c("Numeric"),])>0){

        df_Numeric <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Numeric"),], min)
        df_Numeric <- data.frame(var = DynForest_obj$Inputs$Numeric[df_Numeric$var_split],
                                depth = df_Numeric$depth)

      }else{

        df_Numeric <- data.frame("var" = character(), "depth" = numeric())

      }

      if (nrow(x$V_split[x$V_split$type%in%c("Factor"),])>0){

        df_factor <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Factor"),], min)
        df_factor <- data.frame(var = DynForest_obj$Inputs$Factor[df_factor$var_split],
                                depth = df_factor$depth)

      }else{

        df_factor <- data.frame("var" = character(), "depth" = numeric())

      }

      df <- rbind(df_Numeric, df_factor)

      return(df)

    })

    Other_depth <- aggregate(depth ~ var,
                             do.call(rbind, Other_depth_list), mean)
    Other_depth <- Other_depth[str_order(Other_depth$var),]

    Other_node_depth <- suppressWarnings(
      Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "var",
                                        all.x = TRUE, all.y = TRUE),
             Other_depth_list))
    Other_node_depth <- Other_node_depth[str_order(Other_node_depth$var),]
    colnames(Other_node_depth)[2:ncol(Other_node_depth)] <- paste0("tree", seq(ncol(Other_node_depth)-1))

    # count by tree

    Other_count_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type%in%c("Numeric"),])>0){

        df_Numeric <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Numeric"),], length)
        df_Numeric <- data.frame(var = DynForest_obj$Inputs$Numeric[df_Numeric$var_split],
                                depth = df_Numeric$depth)

      }else{

        df_Numeric <- data.frame("var" = character(), "depth" = numeric())

      }

      if (nrow(x$V_split[x$V_split$type%in%c("Factor"),])>0){

        df_factor <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Factor"),], length)
        df_factor <- data.frame(var = DynForest_obj$Inputs$Factor[df_factor$var_split],
                                depth = df_factor$depth)

      }else{

        df_factor <- data.frame("var" = character(), "depth" = numeric())

      }

      df <- rbind(df_Numeric, df_factor)

      return(df)

    })

    Other_count <- suppressWarnings(
      Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "var",
                                        all.x = TRUE, all.y = TRUE),
             Other_count_list))

    Other_count <- Other_count[str_order(Other_count$var),]
    colnames(Other_count)[2:ncol(Other_count)] <- paste0("tree", seq(ncol(Other_count)-1))

  }else{

    Other_depth <- NULL
    Other_node_depth <- NULL
    Other_count <- NULL

  }

  min_depth <- rbind(Longitudinal_depth, Other_depth)
  min_depth$rank <- rank(min_depth$depth)

  var_node_depth <- rbind(Longitudinal_node_depth, Other_node_depth)
  rownames(var_node_depth) <- seq(nrow(var_node_depth))

  var_count <- rbind(Longitudinal_count, Other_count)
  var_count[is.na(var_count)] <- 0
  rownames(var_count) <- seq(nrow(var_count))

  out <- list(min_depth = min_depth, var_node_depth = var_node_depth,
              var_count = var_count)

  class(out) <- "DynForestVarDepth"

  return(out)

}
