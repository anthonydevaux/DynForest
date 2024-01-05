#' Split function to build the two daughter nodes from longitudinal predictors
#'
#' @param X Input data
#' @param Y Outcome data
#' @param timeVar A character indicating the name of time variable
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param init (Optional) Initial values for linear mixed models
#'
#' @import lcmm
#' @importFrom splines ns
#'
#' @keywords internal
var_split_long <- function(X, Y, timeVar = NULL, nsplit_option = "quantile",
                           cause = 1, nodesize = 1, init = NULL){

  X_ncol <- ncol(X$X)
  all_imp_var <- split_var <- vector("list", X_ncol)
  impur_var <- rep(Inf, X_ncol)
  Pure <- FALSE
  model_param <- list()
  threshold_var <- var_sum <- rep(NA, X_ncol)
  conv_issue <- NULL

  for (i in 1:X_ncol){

    colnames_X_i <- colnames(X$X)[i]

    fixed_var <- all.vars(X$model[[i]]$fixed)
    random_var <- all.vars(X$model[[i]]$random)
    model_var <- unique(c(fixed_var,random_var))

    # compute features from mixed model
    data_model <- data.frame(id = as.numeric(X$id), time = X$time, X$X[,, drop = FALSE])
    colnames(data_model)[which(colnames(data_model)=="time")] <- timeVar
    data_model <- data_model[,c("id", model_var)]

    if (is.null(init[[colnames_X_i]][[1]])){
      init[[colnames_X_i]][[1]] <- NA
    }

    # Mixed model with initial values for parameters ?
    if (!is.na(init[[colnames_X_i]][[1]])){

      model_output <- tryCatch(
        hlme(fixed = X$model[[i]]$fixed,
             random = X$model[[i]]$random,
             subject = "id", data = data_model,
             B = init[[colnames_X_i]],
             maxiter = 100,
             verbose = FALSE),
        error = function(e){ return(NULL) })

      if (is.null(model_output)){ # can occurred with Cholesky matrix inversion

        model_output <- tryCatch(hlme(fixed = X$model[[i]]$fixed,
                                      random = X$model[[i]]$random,
                                      subject = "id", data = data_model,
                                      maxiter = 100,
                                      verbose = FALSE),
                                 error = function(e){ return(NULL) })

      }

    }else{

      model_output <- tryCatch(hlme(fixed = X$model[[i]]$fixed,
                                    random = X$model[[i]]$random,
                                    subject = "id", data = data_model,
                                    maxiter = 100,
                                    verbose = FALSE),
                               error = function(e){ return(NULL) })

    }

    if (model_output$gconv[1]>1e-04 | model_output$gconv[2]>1e-04 | is.null(model_output)){ # convergence issue
      conv_issue <- c(conv_issue, colnames_X_i)
      next()
    }

    init[[colnames_X_i]] <- model_output$best

    model_param[[i]] <- list(beta = model_output$best[(model_output$N[1]+1):model_output$N[2]],
                             varcov = model_output$best[(model_output$N[1]+model_output$N[2]+1):
                                                          (model_output$N[1]+model_output$N[2]+model_output$N[3])],
                             stderr = tail(model_output$best, n = 1),
                             idea0 = model_output$idea0)

    # Random-effect dataframe with NA for subjects where RE cannot be computed
    RE <- merge(unique(Y$id), model_output$predRE, all.x = T, by.x = "x", by.y = "id")[,-1]

    ###########################

    # Other features to compute?

    ###########################

    data_summaries <- RE

    mtry2 <- ncol(data_summaries) # How many features do we draw?
    #mtry_sum <- sample(1:ncol(data_summaries), mtry2)
    mtry_sum <- seq(ncol(data_summaries))

    impur_sum <- rep(Inf, mtry2)
    all_imp_sum <- lapply(seq(mtry2), FUN = function(x) list(Inf, Inf))
    split_sum <- vector("list", mtry2)
    split_sum_threholds <- rep(NA, mtry2)

    for (i_sum in mtry_sum){

      if (!all(is.na(data_summaries[,i_sum]))){

        nsplit <- ifelse(length(unique(na.omit(data_summaries[,i_sum])))>10,
                         10, length(unique(na.omit(data_summaries[,i_sum]))))

        if (nsplit>1) {

          if (nsplit>2){
            if (nsplit_option == "quantile"){
              split_threholds <- unique(quantile(data_summaries[,i_sum], probs = seq(0,1,1/nsplit),
                                                 na.rm = T)[-c(1,nsplit+1)])
            }
            if (nsplit_option == "sample"){
              split_threholds <- unique(sample(data_summaries[,i_sum], nsplit))
            }
          }else{
            split_threholds <- mean(unique(data_summaries[,i_sum]))
          }

          # remove partition according to nodesize criteria
          group_length <- lapply(split_threholds, FUN = function(x){
            table(data_summaries[,i_sum]<=x)
          })

          split_nodesize_ok <- unlist(lapply(group_length, FUN = function(x) !any(x<nodesize)))
          split_threholds <- split_threholds[split_nodesize_ok]
          split_threholds_length <- length(split_threholds)

          if (split_threholds_length>0){ # could happened with tie values

            # Find best feature partition
            split_sum_list <- lapply(seq(split_threholds_length), FUN = function(x){

              split <- ifelse(data_summaries[,i_sum]<=split_threholds[x],1,2)

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

            partition_sum_impur <- unlist(lapply(split_sum_list, function(x) return(x$impur)))

            if (any(partition_sum_impur!=Inf)){
              best_part_sum_nsplit <- which.min(partition_sum_impur)
              split_sum[[i_sum]] <- split_sum_list[[best_part_sum_nsplit]]$split
              impur_sum[i_sum] <- split_sum_list[[best_part_sum_nsplit]]$impur
              all_imp_sum[[i_sum]] <- split_sum_list[[best_part_sum_nsplit]]$imp_list
              split_sum_threholds[i_sum] <- split_threholds[best_part_sum_nsplit]
            }
          }
        }
      }
    }

    if (any(impur_sum!=Inf)){
      best_part_sum <- which.min(impur_sum)
      var_sum[i] <- mtry_sum[best_part_sum]
      split_var[[i]] <- split_sum[[best_part_sum]]
      impur_var[i] <- impur_sum[best_part_sum]
      all_imp_var[[i]] <- all_imp_sum[[best_part_sum]]
      threshold_var[i] <- split_sum_threholds[best_part_sum]
    }

  }

  if (all(impur_var==Inf)){
    return(list(Pure=TRUE))
  }

  var_split <- which.min(impur_var)

  return(list(split = split_var[[var_split]], impur = min(impur_var), impur_list = all_imp_var[[var_split]],
              variable = var_split, variable_summary = var_sum[var_split], threshold = threshold_var[var_split],
              model_param = list(model_param[[var_split]]), conv_issue = conv_issue,
              init = init, Pure = Pure))
}
