#' Mixed-Model Split function
#'
#' @param X
#' @param Y
#' @param timeScale
#'
#' @import kmlShape
#' @import Evomorph
#' @import RiemBase
#' @import lcmm
#' @importFrom splines ns
#'
#' @keywords internal
var_split_MM <- function(X ,Y,timeScale=0.1, nsplit_option = NULL,
                         nodesize = 5){
  # Pour le moment on se concentre sur le cas des variables courbes ::
  impur <- rep(0,dim(X$X)[length(dim(X$X))])
  toutes_imp <- list()
  split <- list()
  centers <- list() # On va stocker les centres associés aux kmeans
  Pure <- FALSE
  model_param <- list()
  threshold <- variable_summary <- rep(NA, ncol(X$X))

  for (i in 1:dim(X$X)[length(dim(X$X))]){

    if (X$type=="factor"){
      if (length(unique(X$X[,i]))>1){
        L <- Fact.partitions(X$X[,i],X$id)
        split_courant <- list()
        impur_courant <- rep(NA,length(L))
        toutes_imp_courant <- list()
        # Il faut maintenant regarder quelles sont les meilleures combinaisons ::
        for (k in 1:length(L)){
          split_courant[[k]] <- rep(2,length(X$id))
          for (l in L[[k]]){
            split_courant[[k]][which(X$id==l)] <- 1
          }
          # Il faut maintenant regarder la qualité du découpage ::
          impurete <- impurity_split(Y,split_courant[[k]])
          impur_courant[k] <- impurete$impur
          toutes_imp_courant[[k]] <- impurete$imp_list
        }
        select <- which.min(impur_courant)
        split[[i]] <- split_courant[[select]]
        impur[i] <- impur_courant[select]
        toutes_imp[[i]] <- toutes_imp_courant[[select]]
      }
      else {
        impur[i] <- Inf
        split[[i]] <- Inf
      }
    }

    if (X$type=="curve"){

      # modele mixte pour construire les resumes

      data_model <- data.frame(id = as.numeric(X$id), time = X$time, marker = X$X[,i])
      colnames(data_model)[which(colnames(data_model)=="marker")] <- colnames(X$X)[i]

      model_output <- hlme(fixed = X$model[[i]]$fixed,
                           random = X$model[[i]]$random,
                           subject = "id", data = data_model, verbose = FALSE)

      model_param[[i]] <- list(beta = model_output$best[(model_output$N[1]+1):model_output$N[2]],
                               varcov = model_output$best[(model_output$N[1]+model_output$N[2]+1):
                                                            (model_output$N[1]+model_output$N[2]+model_output$N[3])],
                               stderr = tail(model_output$best, n = 1),
                               idea0 = model_output$idea0)

      RE <- predRE(model_param[[i]], X$model[[i]], data_model)$bi

      ###########################

      # ici on peut envisager d'autres resumes a calculer

      ###########################

      data_summaries <- RE # on merge tous les resumes

      nsplit <- 10

      mtry2 <- ncol(data_summaries) # nombre de resumes qu'on tire pour chaque variable
      #var_mtry2 <- sample(1:ncol(data_summaries), mtry2)
      var_mtry2 <- seq(ncol(data_summaries))
      
      impurete_sum <- rep(NA, length(var_mtry2))
      split_sum <- list()
      split_threholds_sum <- rep(NA, length(var_mtry2))

      for (i_sum in var_mtry2){ # boucle sur les resumes tires

        if (nsplit_option == "quantile"){ # nsplit sur les quantiles (hors min/max)
          split_threholds <- quantile(data_summaries[,i_sum], probs = seq(0,1,1/nsplit))[-c(1,nsplit+1)]
        }

        if (nsplit_option == "sample"){ # nsplit sur tirage aleatoire d'obversations
          split_threholds <- sample(data_summaries[,i_sum], nsplit)
        }

        impurete_nsplit <- rep(NA, length(split_threholds))
        split_nsplit <- list()

        for (j in 1:length(split_threholds)){ # boucle sur les nsplit

          split_nsplit[[j]] <- factor(ifelse(data_summaries[,i_sum]<=split_threholds[j],1,2))
          impurete <- impurity_split(Y,split_nsplit[[j]], timeScale)
          impurete_nsplit[j] <- impurete$impur

        }

        impurete_sum[i_sum] <- impurete_nsplit[which.min(impurete_nsplit)]
        split_sum[[i_sum]] <- split_nsplit[[which.min(impurete_nsplit)]]
        split_threholds_sum[i_sum] <- split_threholds[which.min(impurete_nsplit)]
      }

      variable_summary[i] <- var_mtry2[which.min(impurete_sum)]
      split[[i]] <- split_sum[[which.min(impurete_sum)]]
      impur[i] <- impurete_sum[which.min(impurete_sum)]
      threshold[i] <- split_threholds_sum[which.min(impurete_sum)]
      toutes_imp[[i]] <- impurete$imp_list # NULL pour surv

    }

    if(X$type=="scalar"){
      if (length(unique(X$X[,i]))>2){

        nsplit <- 10

        if (nsplit_option == "quantile"){ # nsplit sur les quantiles (hors min/max)
          split_threholds <- quantile(X$X[,i], probs = seq(0,1,1/nsplit))[-c(1,nsplit+1)]
        }

        if (nsplit_option == "sample"){ # nsplit sur tirage aleatoire d'obversations
          split_threholds <- sample(X$X[,i], nsplit)
        }

        impurete_nsplit <- rep(NA, length(split_threholds))
        split_nsplit <- list()

        for (j in 1:length(split_threholds)){ # boucle sur les nsplit

          split_nsplit[[j]] <- factor(ifelse(X$X[,i]<=split_threholds[j],1,2))
          impurete <- impurity_split(Y,split_nsplit[[j]], timeScale)
          impurete_nsplit[j] <- impurete$impur

        }

        split[[i]] <- split_nsplit[[which.min(impurete_nsplit)]]
        impur[i] <- impurete_nsplit[which.min(impurete_nsplit)]
        threshold[i] <- split_threholds[which.min(impurete_nsplit)]
        toutes_imp[[i]] <- impurete$imp_list # NULL pour surv

      }

      if (length(unique(X$X[,i]))==2){
        split[[i]] <- rep(2,length(X$X[,i]))
        split[[i]][which(X$X[,i]==unique(X$X[,i])[1])] <- 1
        impurete <- impurity_split(Y,split[[i]], timeScale)
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

  if (Y$type=="surv"){
    if (any(table(split[Y$Y==1],Y$Y[Y$Y==1]) < nodesize)) return(list(Pure=TRUE))
  }

  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur),
              variable_summary=ifelse(X$type=="curve", variable_summary[true_split], NA),
              threshold=ifelse(X$type=="curve"|X$type=="scalar", threshold[true_split], NA),
              model_param=ifelse(X$type=="curve", list(model_param[[true_split]]), NA),
              Pure=Pure))
}
