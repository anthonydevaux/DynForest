#' Internal checking function
#'
#' @param DynForest_obj A \code{DynForest} object resulting from \code{DynForest()} function
#' @param timeData A data.frame containing the id and time measurements variables and the time-dependent predictors.
#' @param fixedData A data.frame containing the id variable and the time-fixed predictors. Non-continuous variables should be characterized as factor.
#' @param idVar A character indicating the name of variable to identify the subjects
#' @param timeVar A character indicating the name of time variable
#' @param timeVarModel A list for each time-dependent predictors containing a list of formula for fixed and random part from the mixed model
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param ntree Number of trees to grow. Default value set to 200.
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#'
#' @keywords internal
checking <- function(DynForest_obj = NULL,
                     timeData, fixedData,
                     idVar, timeVar, timeVarModel,
                     Y, ntree = 200, mtry = 1, nodesize = 1, minsplit = 2){

  if (!inherits(idVar, "character")){
    stop("'idVar' should be a character object!")
  }

  if (!is.null(timeVar)){
    if (!inherits(timeVar, "character")){
      stop("'timeVar' should be a character object!")
    }
  }

  if (is.null(DynForest_obj)){
    if (any(is.null(c(ntree, mtry, nodesize, minsplit)))){
      stop("'ntree', 'mtry', 'nodesize' or 'minsplit' cannot be NULL!")
    }
  }

  if (!is.null(timeData)){
    if (!inherits(timeData, "data.frame")){
      stop("'timeData' should be a data.frame object!")
    }

    if (!any(colnames(timeData)==idVar)){
      stop("'idVar' variable should be contained in timeData!")
    }

    if (!is.null(timeVar)){
      if (!any(colnames(timeData)==timeVar)){
        stop("'timeVar' variable should be contained in timeData!")
      }
    }

    if (is.null(DynForest_obj)){
      if (!inherits(timeVarModel, "list")){
        stop("'timeVarModel' should be a list object!")
      }

      if (!all(colnames(timeData%in%names(timeVarModel)))){
        stop("'timeVarModel' should contain the fixed and random formula for all time-dependent covariates in timeData!")
      }
    }else{

      DynVar <- unlist(DynForest_obj$Inputs)

      if (!all(DynVar%in%c(colnames(timeData),colnames(fixedData)))){
        stop("All variables in DynForest_obj$Inputs should be included in 'timeData' or 'fixedData'!")
      }

    }


    if (!inherits(timeData[,timeVar], c("numeric","integer"))){
      stop("Only continuous time-dependent predictors are allowed in timeData!")
    }

    if (!all(sapply(subset(timeData, select = -c(get(idVar), get(timeVar))),
                    class)%in%c("numeric","integer"))){
      stop("Only continuous time-dependent predictors are allowed in timeData!")
    }

  }

  if (!is.null(fixedData)){
    if (!inherits(fixedData, "data.frame")){
      stop("'fixedData' should be a data.frame object!")
    }

    if (!any(colnames(fixedData)==idVar)){
      stop("'idVar' variable should be contained in fixedData!")
    }

  }

  if (is.null(DynForest_obj)){
    if (is.null(Y)){
      stop("'Y' is missing!")
    }else{
      if (!inherits(Y, "list")){
        stop("'Y' should be a list object!")
      }
      if (!all(names(Y%in%c("type","Y")))){
        stop("'Y' should be a list with type and Y elements!")
      }
      if (!any(colnames(Y$Y)==idVar)){
        stop("'idVar' variable should be contained in Y!")
      }
    }
  }

  if (Y$type=="surv"){
    if (!inherits(Y$Y[,3], "numeric")){
      stop("The column in 'Y$Y' to provide the causes should be typed as numeric with 0 indicating no event!")
    }
  }

}

