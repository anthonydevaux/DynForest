#' Internal checking function
#'
#' @param dynforest_obj A \code{dynforest} object resulting from \code{dynforest()} function
#' @param timeData A data.frame containing the id and time measurements variables and the time-dependent predictors.
#' @param fixedData A data.frame containing the id variable and the time-fixed predictors. Non-continuous variables should be characterized as factor.
#' @param idVar A character indicating the name of variable to identify the subjects
#' @param timeVar A character indicating the name of time variable
#' @param timeVarModel A list for each time-dependent predictors containing a list of formula for fixed and random part from the mixed model
#' @param Y A list of output which should contain: \code{type} defines the nature of the output, can be "\code{surv}", "\code{scalar}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param ntree Number of trees to grow. Default value set to 200.
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#' @param cause (Only with competing events) Number indicates the event of interest.
#'
#' @keywords internal
#' @noRd
checking <- function(dynforest_obj = NULL,
                     timeData, fixedData,
                     idVar, timeVar, timeVarModel,
                     Y, ntree = 200, mtry = 1, nodesize = 1, minsplit = 2,
                     cause = 1){

  # global argument checking
  if (!inherits(idVar, "character")){
    stop("'idVar' should be a character object!")
  }

  if (!is.null(timeVar)){
    if (!inherits(timeVar, "character")){
      stop("'timeVar' should be a character object!")
    }
  }

  if (is.null(dynforest_obj)){
    if (any(is.null(c(ntree, mtry, nodesize, minsplit)))){
      stop("'ntree', 'mtry', 'nodesize' or 'minsplit' cannot be NULL!")
    }
  }

  # timeData checking
  if (!is.null(timeData)){
    if (!inherits(timeData, "data.frame")){
      stop("'timeData' should be a data.frame object!")
    }

    if (!any(colnames(timeData)==idVar)){
      stop("'idVar' variable should be contained in 'timeData'!")
    }

    if (!is.null(timeVar)){
      if (!any(colnames(timeData)==timeVar)){
        stop("'timeVar' variable should be contained in 'timeData'!")
      }
    }

    if (is.null(dynforest_obj)){
      if (!inherits(timeVarModel, "list")){
        stop("'timeVarModel' should be a list object!")
      }

      if (!all(colnames(timeData)[-which(colnames(timeData)%in%c(idVar, timeVar))]%in%names(timeVarModel))){
        stop("'timeData' predictor names should be included in the list names of 'timeVarModel'!")
      }
    }else{

      DynVar <- unlist(dynforest_obj$Inputs)

      if (!all(DynVar%in%c(colnames(timeData),colnames(fixedData)))){
        stop("All variables in dynforest_obj$Inputs should be included in 'timeData' or 'fixedData'!")
      }

    }

    if (!inherits(timeData[,idVar], c("numeric","integer"))){
      stop(paste0(idVar, " variable should be a numeric or integer object in 'timeData'!"))
    }

    if (!inherits(timeData[,timeVar], c("numeric","integer"))){
      stop(paste0(timeVar, " variable should be a numeric or integer object in 'timeData'!"))
    }

    if (!all(sapply(subset(timeData, select = -c(get(idVar), get(timeVar))),
                    class)%in%c("numeric","integer"))){
      stop("Only continuous time-dependent predictors are allowed in 'timeData'!")
    }

  }

  # fixedData checking
  if (!is.null(fixedData)){
    if (!inherits(fixedData, "data.frame")){
      stop("'fixedData' should be a data.frame object!")
    }

    if (!any(colnames(fixedData)==idVar)){
      stop("'idVar' variable should be contained in 'fixedData'!")
    }

    if (!inherits(fixedData[,idVar], c("numeric","integer"))){
      stop(paste0(idVar, " variable should be a numeric or integer object in 'fixedData'!"))
    }

    if (any(duplicated(fixedData[,idVar]))){
      stop("Multiple rows have been found for same id in 'fixedData'!")
    }

  }

  # Y checking
  if (is.null(dynforest_obj)){
    if (is.null(Y)){
      stop("'Y' is missing!")
    }else{
      if (!inherits(Y, "list")){
        stop("'Y' should be a 'list' object!")
      }
      if (!all(names(Y%in%c("type","Y")))){
        stop("'Y' should be a list with 'type' and 'Y' elements!")
      }
      if (!any(colnames(Y$Y)==idVar)){
        stop("'idVar' variable should be contained in Y!")
      }
      if (!inherits(Y$Y[,idVar], c("numeric","integer"))){
        stop(paste0(idVar, " variable should be a numeric or integer object in 'Y$Y'!"))
      }
      if (any(duplicated(Y$Y[,idVar]))){
        stop("Multiple rows have been found for same id in 'Y$Y'!")
      }
      if (Y$type=="surv"){
        if (!inherits(Y$Y[,3], c("numeric","integer"))){
          stop("The column in 'Y$Y' to provide the causes should be typed as numeric with 0 indicating no event!")
        }
        if (all(unique(Y$Y[,3])!=cause)){
          stop("'cause' identifier is not included in 'Y$Y' event column!")
        }
      }
      if (!is.null(fixedData)){
        if (length(fixedData[,idVar])!=length(unique(Y$Y[,idVar]))){
          stop("'fixedData' and 'Y$Y' should contain the same subjects!")
        }
      }
      if (!is.null(timeData)){
        if (length(unique(timeData[,idVar]))!=length(unique(Y$Y[,idVar]))){
          stop("'timeData' and 'Y$Y' should contain the same subjects!")
        }
      }
    }
  }

  # mtry checking
  mtry_max <- ifelse(!is.null(timeData), ncol(timeData)-2, 0) +
    ifelse(!is.null(fixedData), ncol(fixedData)-1, 0)
  if (mtry_max<mtry){
    stop(paste0("'mtry' argument cannot be higher than the maximum allowed (", mtry_max, ")!"))
  }

}

