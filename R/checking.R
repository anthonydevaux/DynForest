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
#' @import cli
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
    cli_abort(c(
      "{.var idVar} must be a character object",
      "x" = "You've supplied a {.cls {class(idVar)}} object"
    ))
  }

  if (!is.null(timeVar)){
    if (!inherits(timeVar, "character")){
      cli_abort(c(
        "{.var timeVar} must be a character object",
        "x" = "You've supplied a {.cls {class(timeVar)}} object"
      ))
    }
  }

  if (is.null(dynforest_obj)){
    if (any(is.null(ntree)|is.null(mtry)|is.null(nodesize)|is.null(minsplit))){
      cli_abort(c(
        "{.var ntree}, {.var mtry}, {.var nodesize} and {.var minsplit} can't be {.var NULL}"
      ))
    }
  }

  # timeData checking
  if (!is.null(timeData)){
    if (!inherits(timeData, "data.frame")){
      cli_abort(c(
        "{.var timeData} must be a data.frame object",
        "x" = "You've supplied a {.cls {class(timeData)}} object"
      ))
    }

    if (!any(colnames(timeData)==idVar)){
      cli_abort(c(
        "Can't find column {.var idVar} in {.var timeData}"
      ))
    }

    if (!is.null(timeVar)){
      if (!any(colnames(timeData)==timeVar)){
        cli_abort(c(
          "Can't find column {.var timeVar} in {.var timeData}"
        ))
      }
    }

    if (is.null(dynforest_obj)){
      if (!inherits(timeVarModel, "list")){
        cli_abort(c(
          "{.var timeVarModel} must be a list object",
          "x" = "You've supplied a {.cls {class(timeVarModel)}} object"
        ))
      }

      if (!all(colnames(timeData)[-which(colnames(timeData)%in%c(idVar, timeVar))]%in%names(timeVarModel))){
        cli_abort(c(
          "{.var timeData} predictor names must be included in the list names of {.var timeVarModel}"
        ))
      }
    }else{

      DynVar <- unlist(dynforest_obj$Inputs)

      if (!all(DynVar%in%c(colnames(timeData),colnames(fixedData)))){
        cli_abort(c(
          "All variables in {.var dynforest_obj$Inputs} must be included in {.var timeData} or {.var fixedData}"
        ))
      }

    }

    if (!inherits(timeData[,idVar], c("numeric","integer"))){
      cli_abort(c(
        "{.var idVar} column in {.var timeData} must be a numeric or integer object",
        "x" = "You've supplied a {.cls {class(timeData[,idVar])}} object"
      ))
    }

    if (!inherits(timeData[,timeVar], c("numeric","integer"))){
      cli_abort(c(
        "{.var timeVar} column in {.var timeData} must be a numeric or integer object",
        "x" = "You've supplied a {.cls {class(timeData[,timeVar])}} object"
      ))
    }

    if (!all(sapply(subset(timeData, select = -c(get(idVar), get(timeVar))),
                    class)%in%c("numeric","integer"))){
      cli_abort(c(
        "Time-dependent predictors in {.var timeVar} must be numeric or integer"
      ))
    }

  }

  # fixedData checking
  if (!is.null(fixedData)){
    if (!inherits(fixedData, "data.frame")){
      cli_abort(c(
        "{.var fixedData} must be a data.frame object",
        "x" = "You've supplied a {.cls {class(fixedData)}} object"
      ))
    }

    if (!any(colnames(fixedData)==idVar)){
      cli_abort(c(
        "Can't find column {.var idVar} in {.var fixedData}"
      ))
    }

    if (!inherits(fixedData[,idVar], c("numeric","integer"))){
      cli_abort(c(
        "{.var idVar} column in {.var fixedData} must be a numeric or integer object",
        "x" = "You've supplied a {.cls {class(fixedData[,idVar])}} object"
      ))
    }

    if (any(duplicated(fixedData[,idVar]))){
      cli_abort(c(
        "Each {.var idVar} identifier in {.var fixedData} must be unique",
      ))
    }

  }

  # Y checking
  if (is.null(dynforest_obj)){
    if (is.null(Y)){
      cli_abort(c(
        "Can't find {.var Y}",
      ))
    }else{
      if (!inherits(Y, "list")){
        cli_abort(c(
          "{.var Y} must be a list object",
          "x" = "You've supplied a {.cls {class(Y)}} object"
        ))
      }
      if (!all(names(Y%in%c("type","Y")))){
        cli_abort(c(
          "{.var Y} must be a list object including {.var type} and {.var Y} elements"
        ))
      }
      if (!any(colnames(Y$Y)==idVar)){
        cli_abort(c(
          "Can't find column {.var idVar} in {.var Y$Y}"
        ))
      }
      if (!inherits(Y$Y[,idVar], c("numeric","integer"))){
        cli_abort(c(
          "{.var idVar} column in {.var Y$Y} must be a numeric or integer object",
          "x" = "You've supplied a {.cls {class(Y$Y[,idVar])}} object"
        ))
      }
      if (any(duplicated(Y$Y[,idVar]))){
        cli_abort(c(
          "Each {.var idVar} identifier in {.var Y$Y} must be unique",
        ))
      }
      if (Y$type=="surv"){
        if (!inherits(Y$Y[,3], c("numeric","integer"))){
          cli_abort(c(
            "Third column (event) in {.var Y$Y} must be a numeric or integer object with 0 indicating no event",
            "x" = "You've supplied a {.cls {class(Y$Y[,3])}} object"
          ))
        }
        if (all(unique(Y$Y[,3])!=cause)){
          cli_abort(c(
            "Can't find any {.var cause} in third column (event) in {.var Y$Y}"
          ))
        }
      }
      if (!is.null(fixedData)){
        if (length(fixedData[,idVar])!=length(unique(Y$Y[,idVar]))){
          cli_abort(c(
            "{.var fixedData} and {.var Y$Y} must contain the same subjects"
          ))
        }
      }
      if (!is.null(timeData)){
        if (length(unique(timeData[,idVar]))!=length(unique(Y$Y[,idVar]))){
          cli_abort(c(
            "{.var timeData} and {.var Y$Y} must contain the same subjects"
          ))
        }
      }
    }
  }

  # mtry checking
  mtry_max <- ifelse(!is.null(timeData), ncol(timeData)-2, 0) +
    ifelse(!is.null(fixedData), ncol(fixedData)-1, 0)
  if (mtry_max<mtry){
    cli_abort(c(
      "{.var mtry} can't be higher than the maximum allowed ({mtry_max})"
    ))
  }

}

