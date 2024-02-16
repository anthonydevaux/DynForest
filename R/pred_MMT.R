#' Predict the leaf by dropping down the subject in the tree
#'
#' @param tree Tree object resulting from \code{Rtmax_surv} function
#' @param Longitudinal A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Numeric A list of numeric predictors which should contain: \code{X} a dataframe with as many columns as numeric predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param timeVar A character indicating the name of time variable
#'
#' @import stringr
#'
#' @keywords internal
pred.MMT <- function(tree, Longitudinal=NULL, Numeric=NULL, Factor=NULL,
                     timeVar = NULL){

  Inputs <- c(Longitudinal$type, Numeric$type, Factor$type)

  id.pred <- unique(get(Inputs[1])$id)

  pred <- rep(NA,length(id.pred))

  pred <- sapply(seq_along(id.pred), FUN = function(i){

    if (is.element("Longitudinal",Inputs)==TRUE) wLongitudinal <- which(Longitudinal$id==id.pred[i])
    if (is.element("Numeric",Inputs)==TRUE) wNumeric <- which(Numeric$id==id.pred[i])
    if (is.element("Factor",Inputs)==TRUE) wFactor <- which(Factor$id==id.pred[i])

    current_node <- 1

    while (is.element(current_node, tree$leaves)==FALSE){

      X <- get(as.character(tree$V_split[which(tree$V_split[,2]==current_node),1]))
      type <- str_to_lower(as.character(tree$V_split[which(tree$V_split[,2]==current_node),1]))
      var.split <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==current_node),3]))
      var.split.sum <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==current_node),4]))
      threshold <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==current_node),5]))

      meanG <- tree$hist_nodes[[as.numeric(2*current_node)]]
      meanD <- tree$hist_nodes[[as.numeric(2*current_node+1)]]

      if (type=="longitudinal"){

        fixed_var <- all.vars(X$model[[var.split]]$fixed)
        random_var <- all.vars(X$model[[var.split]]$random)
        model_var <- unique(c(fixed_var,random_var))

        data_model <- data.frame(id = as.numeric(X$id[wLongitudinal]), time = X$time[wLongitudinal],
                                 X$X[wLongitudinal, , drop = FALSE])
        colnames(data_model)[which(colnames(data_model)=="time")] <- timeVar
        data_model <- data_model[,c("id",model_var)]

        RE <- predRE(tree$model_param[[as.numeric(current_node)]][[1]],
                     X$model[[var.split]], data_model)$bi

        ######################

        # autres resumes

        #####################

        data_summaries <- RE

        if (is.na(data_summaries[,var.split.sum])){
          current_node <- NA
          break
        }

        if (data_summaries[,var.split.sum] < threshold){
          distG <- 0
          distD <- 1
        }else{
          distG <- 1
          distD <- 0
        }

      }
      if (type=="numeric"){

        if (X$X[wNumeric,var.split] < threshold){
          distG <- 0
          distD <- 1
        }else{
          distG <- 1
          distD <- 0
        }

      }

      if (type=="factor"){
        distG <- -1*(is.element(X$X[wFactor,var.split],meanG))
        distD <- -1*(is.element(X$X[wFactor,var.split],meanD))
      }

      if (distG <= distD) { current_node <- 2*current_node}
      if (distD < distG) {current_node <- 2*current_node +1}


    }

    return(current_node)

  })

  return(pred)

}
