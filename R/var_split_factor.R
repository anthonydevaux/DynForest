#' Split function to build the two daughter nodes from factor predictor
#'
#' @param X Input data
#' @param Y Outcome data
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#'
#' @keywords internal
var_split_factor <- function(X, Y, cause = 1, nodesize = 1){

  X_ncol <- ncol(X$X)
  all_imp_var <- split_var <- vector("list", X_ncol)
  impur_var <- rep(Inf, X_ncol)
  Pure <- FALSE

  for (i in 1:X_ncol){

    if (length(unique(X$X[,i]))>1){

      L <- Fact.partitions(X$X[,i],X$id)

      # Find best partition
      split_list <- lapply(seq_along(L), FUN = function(x){

        split <- rep(2,length(X$id))
        split[which(X$id%in%L[[x]])] <- 1

        if ((length(unique(split))>1)&(all(table(split)>=nodesize))){
          # Evaluate the partition
          impur_res <- impurity_split(Y, split, cause = cause)

          impur <- impur_res$impur
          imp_list <- impur_res$imp_list
        }else{
          impur <- Inf
          imp_list <- list(Inf, Inf)
        }

        return(list(split = split, impur = impur, imp_list = imp_list))

      })

      partition_impur <- unlist(lapply(split_list, function(x) return(x$impur)))

      if (any(partition_impur!=Inf)){
        best_part <- which.min(partition_impur)
        split_var[[i]] <- split_list[[best_part]]$split
        impur_var[i] <- split_list[[best_part]]$impur
        all_imp_var[[i]] <- split_list[[best_part]]$imp_list
      }
    }
  }

  if (all(impur_var==Inf)){
    return(list(Pure=TRUE))
  }

  var_split <- which.min(impur_var)

  return(list(split = split_var[[var_split]], impur = min(impur_var), impur_list = all_imp_var[[var_split]],
              variable = var_split, variable_summary = NA, threshold = NA,
              Pure = Pure))
}
