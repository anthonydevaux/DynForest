#' Split function to build the two daughter nodes from numeric predictors
#'
#' @param X Input data
#' @param Y Outcome data
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#'
#' @keywords internal
#' @noRd
var_split_num <- function(X, Y, nsplit_option = "quantile",
                          cause = 1, nodesize = 1){

  X_ncol <- ncol(X$X)
  all_imp_var <- split_var <- vector("list", X_ncol)
  impur_var <- rep(Inf, X_ncol)
  Pure <- FALSE
  threshold_var <- rep(NA, X_ncol)

  for (i in 1:X_ncol){

    if (!all(is.na(X$X[,i]))){

      nsplit <- ifelse(length(unique(na.omit(X$X[,i])))>10,
                       10, length(unique(na.omit(X$X[,i]))))

      if (nsplit>1) {

        if (nsplit>2){
          if (nsplit_option == "quantile"){
            split_threholds <- unique(quantile(X$X[,i], probs = seq(0,1,1/nsplit),
                                               na.rm = T)[-c(1,nsplit+1)])
          }
          if (nsplit_option == "sample"){
            split_threholds <- unique(sample(X$X[,i], nsplit))
          }
        }else{
          split_threholds <- mean(unique(X$X[,i]))
        }

        # remove partition according to nodesize criteria
        group_length <- lapply(split_threholds, FUN = function(x){
          table(X$X[,i]<=x)
        })

        split_nodesize_ok <- unlist(lapply(group_length, FUN = function(x) !any(x<nodesize)))
        split_threholds <- split_threholds[split_nodesize_ok]
        split_threholds_length <- length(split_threholds)

        if (split_threholds_length>0){ # could happened with tie values

          # Find best partition
          split_list <- lapply(seq(split_threholds_length), FUN = function(x){

            split <- ifelse(X$X[,i]<=split_threholds[x],1,2)

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
            threshold_var[i] <- split_threholds[best_part]
          }

        }
      }
    }
  }

  if (all(impur_var==Inf)){
    return(list(Pure=TRUE))
  }

  var_split <- which.min(impur_var)

  return(list(split = split_var[[var_split]], impur = min(impur_var), impur_list = all_imp_var[[var_split]],
              variable = var_split, variable_summary = NA, threshold = threshold_var[var_split], Pure = Pure))
}
