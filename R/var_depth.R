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
#' \dontrun{
#' data(pbc2)
#'
#' # Define time-independent continuous covariate
#' cont_covar <- list(X = pbc2_surv[,"age", drop = FALSE],
#'                   id = pbc2_surv$id)
#'
#' # Define time-independent non continuous covariates
#' fact_covar <- list(X = pbc2_surv[,c("drug","sex")],
#'                    id = pbc2_surv$id)
#'
#' # Define time-dependent continuous markers
#' cont_traj <- list(X = pbc2_long[,c("serBilir","serChol","albumin","alkaline")],
#'                   id = pbc2_long$id,
#'                   time = pbc2_long$time,
#'                   model = list(serBilir = list(fixed = serBilir ~ time,
#'                                                random = ~ time),
#'                                serChol = list(fixed = serChol ~ time + I(time^2),
#'                                               random = ~ time + I(time^2)),
#'                                albumin = list(fixed = albumin ~ time,
#'                                               random = ~ time),
#'                                alkaline = list(fixed = alkaline ~ time,
#'                                                random = ~ time))
#' )
#'
#' # Define outcome (survival here)
#' Y <- list(type = "surv",
#'           Y = Surv(pbc2_surv$years, factor(pbc2_surv$event)),
#'           id = pbc2_surv$id)
#'
#' # Run DynForest function
#' res_dyn <- DynForest(Curve = cont_traj, Factor = fact_covar, Scalar = cont_covar,
#'                      Y = Y, ntree = 200, imp = TRUE,
#'                      mtry = 4, nodesize = 2, minsplit = 3,
#'                      cause = 2)
#'
#' # Run var_depth function
#' res_varDepth <- var_depth(res_dyn)
#'
#' }
#'
#' @export
var_depth <- function(DynForest_obj){

  Inputs <- names(DynForest_obj$Inputs)[unlist(lapply(DynForest_obj$Inputs, FUN = function(x) return(!is.null(x))))]

  if (any(Inputs%in%c("Curve"))){

    # mindepth and node at depth

    Curve_depth_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type=="Curve",])==0){
        return(data.frame("var" = character(), "depth" = numeric()))
      }

      df <- aggregate(depth ~ type + var_split + var_summary, x$V_split[x$V_split$type=="Curve",], min)
      df <- data.frame(var = paste0(DynForest_obj$Inputs$Curve[df$var_split], ".bi", df$var_summary-1),
                       depth = df$depth)
      return(df)

    })

    Curve_depth <- aggregate(depth ~ var,
                             do.call(rbind, Curve_depth_list), mean)
    Curve_depth <- Curve_depth[str_order(Curve_depth$var),]

    Curve_node_depth <- suppressWarnings(
      Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "var",
                                        all.x = TRUE, all.y = TRUE),
             Curve_depth_list))
    Curve_node_depth <- Curve_node_depth[str_order(Curve_node_depth$var),]
    colnames(Curve_node_depth)[2:ncol(Curve_node_depth)] <- paste0("tree", seq(ncol(Curve_node_depth)-1))

    # count by tree

    Curve_count_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type=="Curve",])==0){
        return(data.frame("var" = character(), "depth" = numeric()))
      }

      df <- aggregate(depth ~ type + var_split + var_summary, x$V_split[x$V_split$type=="Curve",], length)
      df <- data.frame(var = paste0(DynForest_obj$Inputs$Curve[df$var_split], ".bi", df$var_summary-1),
                       depth = df$depth)
      return(df)

    })

    Curve_count <- suppressWarnings(
      Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "var",
                                        all.x = TRUE, all.y = TRUE),
             Curve_count_list))

    Curve_count <- Curve_count[str_order(Curve_count$var),]
    colnames(Curve_count)[2:ncol(Curve_count)] <- paste0("tree", seq(ncol(Curve_count)-1))

  }else{

    Curve_depth <- NULL
    Curve_node_depth <- NULL
    Curve_count <- NULL

  }

  if (any(Inputs%in%c("Scalar","Factor"))){

    # mindepth and node at depth

    Other_depth_list <- apply(DynForest_obj$rf, 2, FUN = function(x){

      if (nrow(x$V_split[x$V_split$type%in%c("Scalar"),])>0){

        df_scalar <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Scalar"),], min)
        df_scalar <- data.frame(var = DynForest_obj$Inputs$Scalar[df_scalar$var_split],
                                depth = df_scalar$depth)

      }else{

        df_scalar <- data.frame("var" = character(), "depth" = numeric())

      }

      if (nrow(x$V_split[x$V_split$type%in%c("Factor"),])>0){

        df_factor <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Factor"),], min)
        df_factor <- data.frame(var = DynForest_obj$Inputs$Factor[df_factor$var_split],
                                depth = df_factor$depth)

      }else{

        df_factor <- data.frame("var" = character(), "depth" = numeric())

      }

      df <- rbind(df_scalar, df_factor)

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

      if (nrow(x$V_split[x$V_split$type%in%c("Scalar"),])>0){

        df_scalar <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Scalar"),], length)
        df_scalar <- data.frame(var = DynForest_obj$Inputs$Scalar[df_scalar$var_split],
                                depth = df_scalar$depth)

      }else{

        df_scalar <- data.frame("var" = character(), "depth" = numeric())

      }

      if (nrow(x$V_split[x$V_split$type%in%c("Factor"),])>0){

        df_factor <- aggregate(depth ~ type + var_split, x$V_split[x$V_split$type%in%c("Factor"),], length)
        df_factor <- data.frame(var = DynForest_obj$Inputs$Factor[df_factor$var_split],
                                depth = df_factor$depth)

      }else{

        df_factor <- data.frame("var" = character(), "depth" = numeric())

      }

      df <- rbind(df_scalar, df_factor)

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

  min_depth <- rbind(Curve_depth, Other_depth)
  min_depth$rank <- rank(min_depth$depth)

  var_node_depth <- rbind(Curve_node_depth, Other_node_depth)
  rownames(var_node_depth) <- seq(nrow(var_node_depth))

  var_count <- rbind(Curve_count, Other_count)
  var_count[is.na(var_count)] <- 0
  rownames(var_count) <- seq(nrow(var_count))

  return(list(min_depth = min_depth, var_node_depth = var_node_depth,
              var_count = var_count))

}
