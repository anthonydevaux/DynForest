#' Grow random survival tree using multivariate longitudinal endogenous covariates
#'
#' @param Y A list of output which should contain: \code{type} defines the nature of the outcome, can be "\code{surv}", "\code{numeric}" or "\code{factor}"; \code{Y} is the output variable; \code{id} is the vector of the identifiers for each individuals, they should be the same as the identifiers of the inputs.
#' @param Longitudinal A list of longitudinal predictors which should contain: \code{X} a dataframe with one row for repeated measurement and as many columns as markers; \code{id} is the vector of the identifiers for the repeated measurements contained in \code{X}; \code{time} is the vector of the measurement times contained in \code{X}.
#' @param Numeric A list of numeric predictors which should contain: \code{X} a dataframe with as many columns as numeric predictors; \code{id} is the vector of the identifiers for each individual.
#' @param Factor A list of factor predictors which should contain: \code{X} a dataframe with as many columns as factor predictors; \code{id} is the vector of the identifiers for each individual.
#' @param timeVar A character indicating the name of time variable
#' @param mtry Number of candidate variables randomly drawn at each node of the trees. This parameter should be tuned by minimizing the OOB error. Default is `NULL`.
#' @param nsplit_option A character indicates how the values are chosen to build the two groups for the splitting rule (only for continuous predictors). Values are chosen using deciles (\code{nsplit_option}="quantile") or randomly (\code{nsplit_option}="sample"). Default value is "quantile".
#' @param nodesize Minimal number of subjects required in both child nodes to split. Cannot be smaller than 1.
#' @param minsplit (Only with survival outcome) Minimal number of events required to split the node. Cannot be smaller than 2.
#' @param cause (Only with competing events) Number indicates the event of interest.
#' @param seed Seed to replicate results
#'
#' @import stringr
#' @import survival
#' @import prodlim
#' @importFrom splines ns
#'
#' @keywords internal
#' @noRd
DynTree_surv <- function(Y, Longitudinal = NULL, Numeric = NULL, Factor = NULL,
                         timeVar = NULL, mtry = 1, nsplit_option = "quantile",
                         nodesize = 1, minsplit = 2, cause = 1, seed = 1234){

  Inputs <- c(Longitudinal$type, Numeric$type, Factor$type)
  type_pred <- unlist(sapply(Inputs, FUN = function(x) return(rep(get(x)$type, ncol(get(x)$X)))))

  V_split <- data.frame(type = character(), id_node = integer(), var_split = integer(),
                        feature = integer(), threshold = numeric(), N = integer(),
                        Nevent = integer(), stringsAsFactors = FALSE)

  # Bootstrap sample
  set.seed(seed)
  id_boot <- unique(sample(unique(Y$id), length(unique(Y$id)), replace=TRUE))

  # Longitudinal bootstrap data
  if (!is.null(Longitudinal)){
    wXLongitudinal <- which(Longitudinal$id%in%id_boot)
    Longitudinal_boot <- list(type = Longitudinal$type,
                              X = Longitudinal$X[wXLongitudinal,, drop=FALSE],
                              id = Longitudinal$id[wXLongitudinal],
                              time = Longitudinal$time[wXLongitudinal],
                              model = Longitudinal$model)
  }

  # Numeric bootstrap data
  if (!is.null(Numeric)){
    wXNumeric <- which(Numeric$id%in%id_boot)
    Numeric_boot <- list(type = Numeric$type,
                         X = Numeric$X[wXNumeric,, drop=FALSE],
                         id = Numeric$id[wXNumeric])
  }

  # Factor bootstrap data
  if (!is.null(Factor)){
    wXFactor <- which(Factor$id%in%id_boot)
    Factor_boot <- list(type = Factor$type,
                        X = Factor$X[wXFactor,, drop=FALSE],
                        id = Factor$id[wXFactor])
  }

  # Outcome bootstrap data
  wY <- which(Y$id%in%id_boot)
  Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY], comp=Y$comp)


  # Initialize the tree
  id_nodes <- rep(1,length(Y_boot$id)) # nodes
  id_leaves <- NULL
  current_nodes <- 1
  Y_pred <- hist_nodes <- list()

  # Initialize mixed models lists
  model_init <- model_param <- conv_issue <- list()

  for (p in seq_along(unique(Y_boot$id)/2-1)){

    for (current_node in current_nodes){

      current_node_chr <- as.character(current_node)

      # mtry predictors
      set.seed(seed+p*which(current_node==current_nodes))
      mtry_pred <- sample(type_pred, mtry)
      mtry_type_pred <- unique(mtry_pred)

      # Id data at current node
      w <- which(id_nodes==current_node)
      unique_Y_boot_id_w <- unique(Y_boot$id[w])

      if (!is.null(Longitudinal)) wXLongitudinal <- which(Longitudinal_boot$id%in%unique_Y_boot_id_w)
      if (!is.null(Numeric)) wXNumeric <- which(Numeric_boot$id%in%unique_Y_boot_id_w)
      if (!is.null(Factor)) wXFactor <- which(Factor_boot$id%in%unique_Y_boot_id_w)

      Y_current <- list(type=Y_boot$type, Y=Y_boot$Y[w], id=Y_boot$id[w], comp=Y$comp)
      Nevent_current <- sum(Y_current$Y[,2]==cause)
      N_current <- length(Y_current$id)

      F_SPLIT <- data.frame(TYPE = character(), Impurity = numeric(), stringsAsFactors = FALSE)
      leaf_flag <- FALSE

      isLongitudinal <- is.element("Longitudinal", mtry_type_pred)
      isNumeric <- is.element("Numeric", mtry_type_pred)
      isFactor <- is.element("Factor", mtry_type_pred)

      # Node can be split?
      if (Nevent_current >= minsplit && N_current >= nodesize*2){

        # Data at current node with mtry predictors
        if (isLongitudinal){

          tirageLongitudinal <- sample(1:ncol(Longitudinal$X),length(which(mtry_pred=="Longitudinal")))
          Longitudinal_current <- list(type = Longitudinal_boot$type,
                                       X=Longitudinal_boot$X[wXLongitudinal,tirageLongitudinal, drop=FALSE],
                                       id=Longitudinal_boot$id[wXLongitudinal, drop=FALSE],
                                       time=Longitudinal_boot$time[wXLongitudinal, drop=FALSE],
                                       model = Longitudinal_boot$model[tirageLongitudinal])

          if (current_node > 1){
            model_init <- getParamMM(current_node = current_node, markers = colnames(Longitudinal_current$X),
                                     params = model_init)
          }else{
            model_init[[current_node_chr]] <- lapply(Longitudinal$model, FUN = function(x) x$init.param)
          }

        }

        if (isNumeric){
          tirageNumeric <- sample(1:ncol(Numeric$X),length(which(mtry_pred=="Numeric")))
          Numeric_current <- list(type = Numeric_boot$type, X=Numeric_boot$X[wXNumeric,tirageNumeric, drop=FALSE],
                                  id=Numeric_boot$id[wXNumeric, drop=FALSE])
        }

        if (isFactor){
          tirageFactor <- sample(1:ncol(Factor$X),length(which(mtry_pred=="Factor")))
          Factor_current <- list(type = Factor_boot$type, X=Factor_boot$X[wXFactor,tirageFactor, drop=FALSE],
                                 id=Factor_boot$id[wXFactor, drop=FALSE])
        }

        # Try best split on mtry predictors
        if (is.element("Factor", mtry_type_pred)){

          leaf_split_Factor <- var_split_factor(X = Factor_current, Y = Y_current,
                                                cause = cause, nodesize = nodesize)

          if (leaf_split_Factor$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,
                             data.frame(TYPE = "Factor", Impurity = leaf_split_Factor$impur,
                                        stringsAsFactors = FALSE))
          }
        }

        if (is.element("Longitudinal", mtry_type_pred)){

          leaf_split_Longitudinal <- var_split_long(X = Longitudinal_current, Y = Y_current,
                                                    timeVar = timeVar,
                                                    nsplit_option = nsplit_option,
                                                    cause = cause, nodesize = nodesize,
                                                    init = model_init[[current_node_chr]])

          if (leaf_split_Longitudinal$Pure==FALSE){
            model_init[[current_node_chr]] <- leaf_split_Longitudinal$init # update initial values at current node
            F_SPLIT <- rbind(F_SPLIT,
                             data.frame(TYPE = "Longitudinal", Impurity = leaf_split_Longitudinal$impur,
                                        stringsAsFactors = FALSE))

            conv_issue[[current_node_chr]] <- leaf_split_Longitudinal$conv_issue
          }
        }

        if (is.element("Numeric", mtry_type_pred)){

          leaf_split_Numeric <- var_split_num(X = Numeric_current, Y = Y_current,
                                              nsplit_option = nsplit_option,
                                              cause = cause, nodesize = nodesize)

          if (leaf_split_Numeric$Pure==FALSE){
            F_SPLIT <- rbind(F_SPLIT,
                             data.frame(TYPE = "Numeric", Impurity = leaf_split_Numeric$impur,
                                        stringsAsFactors = FALSE))
          }


        }

      }else{
        leaf_flag <- TRUE
      }

      if (nrow(F_SPLIT)>0){

        best_split_type <- F_SPLIT$TYPE[which.min(F_SPLIT$Impurity)]
        X_boot <- get(paste0(best_split_type, "_current"))

        # Get best partition
        leaf_split <- get(paste0("leaf_split_", best_split_type))
        best_pred <- get(paste0("tirage", best_split_type))[leaf_split$variable]

        left_id <- unique(Y_current$id)[which(leaf_split$split==1)]
        right_id <- unique(Y_current$id)[which(leaf_split$split==2)]

        length_left <- length(left_id)
        length_right <- length(right_id)

        if (length_left<nodesize | length_right<nodesize){
          leaf_flag <- TRUE
        }

      }else{
        leaf_flag <- TRUE
      }

      if (!leaf_flag){

        # add node split to V_split
        V_split <- rbind(V_split,
                         data.frame(type = best_split_type, id_node = current_node,
                                    var_split = best_pred, feature = leaf_split$variable_summary,
                                    threshold = leaf_split$threshold, N = N_current,
                                    Nevent = Nevent_current, stringsAsFactors = FALSE))

        model_param[[current_node_chr]] <- leaf_split$model_param

        w_left <- which(X_boot$id%in%left_id)
        wY_left <- which(Y_boot$id%in%left_id)

        w_right <- which(X_boot$id%in%right_id)
        wY_right <- which(Y_boot$id%in%right_id)

        # Check for missing split
        if (anyNA(leaf_split$split)){
          na_id <- unique(Y_current$id)[which(is.na(leaf_split$split))]
        }else{
          na_id <- NULL
        }

        if (!is.null(na_id)){
          wY_na <- which(Y_boot$id%in%na_id)
          id_nodes[wY_na] <- NA
        }

        id_nodes[wY_left] <- 2*current_node
        id_nodes[wY_right] <- 2*current_node+1

        if (best_split_type=="Longitudinal"){
          meanFg <- NA
          meanFd <- NA
        }

        if (best_split_type=="Factor"){
          meanFg <- unique(X_boot$X[w_left, leaf_split$variable])
          meanFd <- unique(X_boot$X[w_right, leaf_split$variable])
        }

        if (best_split_type=="Numeric"){
          meanFg <- mean(X_boot$X[w_left, leaf_split$variable])
          meanFd <- mean(X_boot$X[w_right, leaf_split$variable])
        }

        hist_nodes[[as.character(2*current_node)]] <- meanFg
        hist_nodes[[as.character(2*current_node+1)]] <- meanFd

      }else{

        id_leaves <- c(id_leaves, current_node)

        V_split <- rbind(V_split,
                         data.frame(type = "Leaf", id_node = current_node, var_split = NA,
                                    feature = NA, threshold = NA, N = N_current,
                                    Nevent = Nevent_current, stringsAsFactors = FALSE))

      }

    }

    current_nodes <- setdiff(unique(na.omit(id_nodes)), id_leaves)

  }

  # depth level
  V_split <- V_split[order(V_split$id_node),]
  V_split$depth <- floor(log(V_split$id_node, base = 2)) + 1

  # Get prediction for each leaf
  if (nrow(V_split)>0){
    rownames(V_split) <- seq(nrow(V_split))
  }

  for (q in sort(unique(na.omit(id_nodes)))){

    w <- which(id_nodes == q)

    datasurv <- data.frame(time_event = Y_boot$Y[w][,1], event = Y_boot$Y[w][,2])
    fit <- prodlim(Hist(time_event, event)~1, data = datasurv,
                   type = "risk")

    if (is.null(fit$cuminc)){

      if (all(unique(datasurv$event)==0)){ # case with no event
        pred <- list()
        pred[[as.character(cause)]] <- data.frame(times=fit$time, traj = 0) # no event => no risk
      }else{ # case with event no matter which one
        u_current_causes <- unique(datasurv$event)
        current_cause <- u_current_causes[u_current_causes!=0] # keep leaf cause
        pred <- list(data.frame(times=fit$time, traj=1-fit$surv)) # 1-KM
        names(pred) <- as.character(current_cause)
      }

    }else{
      pred <- lapply(fit$cuminc, FUN = function(x) return(data.frame(times=fit$time, traj=x))) # CIF Aalen-Johansen
    }

    Y_pred[[as.character(q)]] <- lapply(pred, function(x){
      combine_times(pred = x, newtimes = unique(Y$Y[,1]), type = "risk")
    })

  }

  return(list(leaves = id_nodes, idY = Y_boot$id, Ytype = Y_boot$type, V_split = V_split,
              hist_nodes = hist_nodes, Y_pred = Y_pred, Y = Y, boot = id_boot, conv_issue = conv_issue,
              model_param = model_param))

}
