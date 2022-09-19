#' Predict the leaf by dropping down the subject in the tree
#'
#' @param tree Tree object resulting from \code{Rtmax_surv} function
#' @param Curve A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Scalar A list of scalar predictors which should contain: \code{X} a dataframe with as many columns as scalar predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#'
#' @import stringr
#'
#' @keywords internal
pred.MMT <- function(tree, Curve=NULL,Scalar=NULL,Factor=NULL){

  Inputs <- read.Xarg(c(Curve,Scalar,Factor))

  id.pred <- unique(get(Inputs[1])$id)

  pred <- rep(NA,length(id.pred))

  for (i in 1:length(id.pred)){

    if (is.element("Curve",Inputs)==TRUE) wCurve <- which(Curve$id==id.pred[i])
    if (is.element("Scalar",Inputs)==TRUE) wScalar <- which(Scalar$id==id.pred[i])
    if (is.element("Factor",Inputs)==TRUE) wFactor <- which(Factor$id==id.pred[i])

    noeud_courant <- 1

    while (is.element(noeud_courant, tree$feuilles)==FALSE){

      X <- get(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      type <- str_to_lower(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),1]))
      var.split <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),3]))
      var.split.sum <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),4]))
      threshold <- as.numeric(as.character(tree$V_split[which(tree$V_split[,2]==noeud_courant),5]))

      meanG <- tree$hist_nodes[[2*noeud_courant]]
      meanD <- tree$hist_nodes[[2*noeud_courant+1]]

      if (type=="curve"){

        fixed_var <- all.vars(X$model[[var.split]]$fixed)
        random_var <- all.vars(X$model[[var.split]]$random)
        model_var <- unique(c(fixed_var,random_var))

        data_model <- data.frame(id = as.numeric(X$id[wCurve]), time = X$time[wCurve],
                                 X$X[wCurve, , drop = FALSE])

        data_model <- data_model[,c("id",model_var)]

        RE <- predRE(tree$model_param[[noeud_courant]][[1]],
                     X$model[[var.split]], data_model)$bi

        ######################

        # autres resumes

        #####################

        data_summaries <- RE

        if (is.na(data_summaries[,var.split.sum])){
          noeud_courant <- NA
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
      if (type=="scalar"){

        if (X$X[wScalar,var.split] < threshold){
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

      if (distG <= distD) { noeud_courant <- 2*noeud_courant}
      if (distD < distG) {noeud_courant <- 2*noeud_courant +1}


    }

    pred[i] <- noeud_courant

  }
  return(pred)
}
