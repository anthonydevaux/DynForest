#' Split function to build the two daughter nodes
#'
#' @param X Input data
#' @param Y Outcome data
#' @param timeVar A character indicating the name of time variable
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param init (Optional) Initial values for linear mixed models
#'
#' @import lcmm
#' @importFrom splines ns
#'
#' @keywords internal
var_split <- function(X, Y, timeVar = NULL, nsplit_option = "quantile",
                      nodesize = 1, init = NULL){

  impur <- rep(0,ncol(X$X))
  toutes_imp <- list()
  split <- list()
  Pure <- FALSE
  model_param <- list()
  threshold <- variable_summary <- rep(NA, ncol(X$X))
  impurete <- NULL

  for (i in 1:ncol(X$X)){

    if (X$type=="Factor"){
      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()
        # Find best partition
        for (k in 1:length(L)){

          split_courant[[k]] <- rep(2,length(X$id))
          split_courant[[k]][which(X$id%in%L[[k]])] <- 1

          if ((length(unique(split_courant[[k]]))==1)|(any(table(split_courant[[k]])<nodesize))){
            impur_courant[k] <- Inf
            next()
          }

          # Evaluate the partition
          impurete <- impurity_split(Y, split_courant[[k]])
          impur_courant[k] <- impurete$impur
          toutes_imp_courant[[k]] <- impurete$imp_list
        }

        if (!is.null(impurete) & !all(impur_courant==Inf)){
          select <- which.min(impur_courant)
          split[[i]] <- split_courant[[select]]
          impur[i] <- impur_courant[select]
          toutes_imp[[i]] <- toutes_imp_courant[[select]]
        }else{
          split[[i]] <- Inf
          impur[i] <- Inf
        }

      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }

    if (X$type=="Longitudinal"){

      fixed_var <- all.vars(X$model[[i]]$fixed)
      random_var <- all.vars(X$model[[i]]$random)
      model_var <- unique(c(fixed_var,random_var))

      # compute features from mixed model
      data_model <- data.frame(id = as.numeric(X$id), time = X$time, X$X[,, drop = FALSE])
      colnames(data_model)[which(colnames(data_model)=="time")] <- timeVar
      data_model <- data_model[,c("id", model_var)]

      if (is.null(init[[colnames(X$X)[i]]][[1]])){
        init[[colnames(X$X)[i]]][[1]] <- NA
      }

      # Mixed model with initial values for parameters ?
      if (!is.na(init[[colnames(X$X)[i]]][[1]])){

        model_output <- tryCatch(
          hlme(fixed = X$model[[i]]$fixed,
               random = X$model[[i]]$random,
               subject = "id", data = data_model,
               B = init[[colnames(X$X)[i]]],
               maxiter = 100,
               verbose = FALSE),
          error = function(e){ return(NULL) })

        if (is.null(model_output)){ # can occurred with Cholesky matrix inversion

          model_output <- hlme(fixed = X$model[[i]]$fixed,
                               random = X$model[[i]]$random,
                               subject = "id", data = data_model,
                               maxiter = 100,
                               verbose = FALSE)

        }

      }else{

        model_output <- hlme(fixed = X$model[[i]]$fixed,
                             random = X$model[[i]]$random,
                             subject = "id", data = data_model,
                             maxiter = 100,
                             verbose = FALSE)

      }

      if (model_output$gconv[1]>1e-04 | model_output$gconv[2]>1e-04){ # convergence issue

        impur[i] <- Inf
        split[[i]] <- Inf
        next()

      }

      init[[colnames(X$X)[i]]] <- model_output$best

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
      #var_mtry2 <- sample(1:ncol(data_summaries), mtry2)
      var_mtry2 <- seq(ncol(data_summaries))

      impurete_sum <- rep(NA, length(var_mtry2))
      split_sum <- impurete_imp_list_sum <- vector("list", length = length(var_mtry2))
      split_threholds_sum <- rep(NA, length(var_mtry2))

      for (i_sum in var_mtry2){

        if (!all(is.na(data_summaries[,i_sum]))){

          nsplit <- ifelse(length(unique(na.omit(data_summaries[,i_sum])))>10,
                                  10, length(unique(na.omit(data_summaries[,i_sum]))))

          if (nsplit==1) {
            impurete_sum[i_sum] <- NA
            split_sum[[i_sum]] <- NA
            impurete_imp_list_sum[[i_sum]] <- NA
            split_threholds_sum[i_sum] <- NA
            next()
          }

          if (nsplit_option == "quantile"){
            split_threholds <- unique(quantile(data_summaries[,i_sum], probs = seq(0,1,1/nsplit),
                                               na.rm = T)[-c(1,nsplit+1)])
          }

          if (nsplit_option == "sample"){
            split_threholds <- unique(sample(data_summaries[,i_sum], nsplit))
          }

          # remove partition according to nodesize criteria
          group_length <- lapply(split_threholds, FUN = function(x){
            table(data_summaries[,i_sum]<=x)
          })

          split_nodesize_ok <- unlist(lapply(group_length, FUN = function(x) !any(x<nodesize)))
          split_threholds <- split_threholds[split_nodesize_ok]

          if (length(split_threholds)==0){ # could happened with tie values
            impurete_sum[i_sum] <- NA
            split_sum[[i_sum]] <- NA
            impurete_imp_list_sum[[i_sum]] <- NA
            split_threholds_sum[i_sum] <- NA
            next()
          }

          impurete_nsplit <- rep(NA, length(split_threholds))
          split_nsplit <- impurete_imp_list_nsplit <- vector("list", length = length(split_threholds))

          for (j in 1:length(split_threholds)){

            split_nsplit[[j]] <- factor(ifelse(data_summaries[,i_sum]<=split_threholds[j],1,2))

            if (length(unique(split_nsplit[[j]]))==1){
              impurete_nsplit[j] <- Inf
              next()
            }

            impurete <- impurity_split(Y, split_nsplit[[j]])
            impurete_nsplit[j] <- impurete$impur
            impurete_imp_list_nsplit[[j]] <- impurete$imp_list

          }

          if (!is.null(impurete) & !all(impurete_nsplit==Inf)){
            impurete_sum[i_sum] <- impurete_nsplit[which.min(impurete_nsplit)]
            split_sum[[i_sum]] <- split_nsplit[[which.min(impurete_nsplit)]]
            impurete_imp_list_sum[[i_sum]] <- impurete_imp_list_nsplit[[which.min(impurete_nsplit)]]
            split_threholds_sum[i_sum] <- split_threholds[which.min(impurete_nsplit)]
          }else{
            impurete_sum[i_sum] <- NA
            split_sum[[i_sum]] <- NA
            impurete_imp_list_sum[[i_sum]] <- NA
            split_threholds_sum[i_sum] <- NA
          }

        }else{

          # if random effects cant be computed
          impurete_sum[i_sum] <- NA
          split_sum[[i_sum]] <- NA
          impurete_imp_list_sum[[i_sum]] <- NA
          split_threholds_sum[i_sum] <- NA

        }

      }

      if (length(which.min(impurete_sum))>0){
        variable_summary[i] <- var_mtry2[which.min(impurete_sum)]
        split[[i]] <- split_sum[[which.min(impurete_sum)]]
        impur[i] <- impurete_sum[which.min(impurete_sum)]
        threshold[i] <- split_threholds_sum[which.min(impurete_sum)]
        toutes_imp[[i]] <- impurete_imp_list_sum[[which.min(impurete_sum)]]
      }else{
        impur[i] <- Inf
        split[[i]] <- Inf
      }

    }

    if (X$type=="Numeric"){
      if (length(unique(X$X[,i]))>2){

        nsplit <- ifelse(length(unique(X$X[,i]))>10, 10, length(unique(X$X[,i])))

        if (nsplit_option == "quantile"){
          split_threholds <- unique(quantile(X$X[,i], probs = seq(0,1,1/nsplit),
                                             na.rm = T)[-c(1,nsplit+1)])
        }

        if (nsplit_option == "sample"){
          split_threholds <- unique(sample(X$X[,i], nsplit))
        }

        # remove partition according to nodesize criteria
        group_length <- lapply(split_threholds, FUN = function(x){
          table(X$X[,i]<=x)
        })

        split_nodesize_ok <- unlist(lapply(group_length, FUN = function(x) !any(x<nodesize)))
        split_threholds <- split_threholds[split_nodesize_ok]

        if (length(split_threholds)>0){ # could happened with tie values

          impurete_nsplit <- rep(NA, length(split_threholds))
          split_nsplit <- impurete_imp_list_nsplit <- vector("list", length = length(split_threholds))

          for (j in 1:length(split_threholds)){

            split_nsplit[[j]] <- factor(ifelse(X$X[,i]<=split_threholds[j],1,2))

            if (length(unique(split_nsplit[[j]]))==1){
              impurete_nsplit[j] <- Inf
              next()
            }

            impurete <- impurity_split(Y, split_nsplit[[j]])
            impurete_nsplit[j] <- impurete$impur
            impurete_imp_list_nsplit[[j]] <- impurete$imp_list

          }
        }else{
          impurete_nsplit <- Inf
        }

        if (!is.null(impurete) & !all(impurete_nsplit==Inf)){
          split[[i]] <- split_nsplit[[which.min(impurete_nsplit)]]
          impur[i] <- impurete_nsplit[which.min(impurete_nsplit)]
          threshold[i] <- split_threholds[which.min(impurete_nsplit)]
          toutes_imp[[i]] <- impurete_imp_list_nsplit[[which.min(impurete_nsplit)]]
        }else{
          impur[i] <- Inf
          split[[i]] <- Inf
        }

      }

      if (length(unique(X$X[,i]))==2){
        split[[i]] <- rep(2,length(X$X[,i]))
        split[[i]][which(X$X[,i]==unique(X$X[,i])[1])] <- 1
        impurete <- impurity_split(Y, split[[i]])
        impur[i] <- impurete$impur
        threshold[i] <- mean(unique(X$X[,i]))
        toutes_imp[[i]] <- impurete$imp_list
      }

      if (length(unique(X$X[,i]))==1) {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }
  }

  if (length(unique(impur))==1 & is.element(Inf,impur)==TRUE){
    return(list(Pure=TRUE))
  }
  true_split <- which.min(impur)
  split <- split[[true_split]]

  return(list(split=split, impur=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur),
              variable_summary=ifelse(X$type=="Longitudinal", variable_summary[true_split], NA),
              threshold=ifelse(X$type=="Longitudinal"|X$type=="Numeric", threshold[true_split], NA),
              model_param=ifelse(X$type=="Longitudinal", list(model_param[[true_split]]), NA),
              init = init,
              Pure=Pure))
}
