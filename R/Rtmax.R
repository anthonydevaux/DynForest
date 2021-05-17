#' Randomized Dynamic tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Y [list]:
#' @param mtry [integer]:
#' @param timeScale [numeric]:
#' @param nsplit_option
#' @param nodesize
#' @param cause
#'
#' @import kmlShape
#' @import RiemBase
#' @import stringr
#' @import Evomorph
#' @import geomorph
#' @import survival
#' @import prodlim
#' @importFrom splines ns
#'
#' @keywords internal
Rtmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Y=NULL, mtry = 1, timeScale = 0.1,
                  nsplit_option = "quantile", nodesize = 1, cause = 1){

  inputs <- read.Xarg(c(Curve,Scalar,Factor))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  impurity_feuilles <- NULL
  V_split <- data.frame(type = character(), num_noeud = integer(), var_split = integer(),
                        var_summary = integer(), threshold = numeric(), N = integer(),
                        Nevent = integer(), stringsAsFactors = FALSE)
  hist_nodes <- list()
  model_param <- list()
  id_boot <- unique(sample(unique(Y$id), length(unique(Y$id)), replace=TRUE))
  boot <- id_boot
  decoupe <- 1

  wXCurve <- NULL
  wXScalar <- NULL
  wXFactor <- NULL
  wY <- NULL

  for (k in id_boot){
    wY <- c(wY, which(Y$id==k))
    if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve$id==k))
    if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar$id==k))
    if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor$id==k))
  }

  Y_pred <- list()

  if (is.element("curve",inputs)==TRUE) Curve_boot <- list(type=Curve$type,   X=Curve$X[wXCurve,, drop=FALSE], id= Curve$id[wXCurve], time = Curve$time[wXCurve],
                                                           model=Curve$model) ### bootstrap pour les courbes
  if (is.element("scalar",inputs)==TRUE) Scalar_boot <- list(type=Scalar$type,   X=Scalar$X[wXScalar,, drop=FALSE], id= Scalar$id[wXScalar]) ### bootstrap pour les courbes
  if (is.element("factor",inputs)==TRUE) Factor_boot <- list(type=Factor$type,   X=Factor$X[wXFactor,, drop=FALSE], id= Factor$id[wXFactor])


  if (Y$type=="curve") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY], time=Y$time[wY])} ### idem pour Y
  if (Y$type=="surv") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY], comp=Y$comp)}
  if (Y$type=="factor" || Y$type=="scalar") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY])}


  imp_nodes <- list()
  imp_nodes[[1]] = Inf
  impurete = Inf
  if (Y$type!="surv"){
    impurete <- impurity(Y_boot, timeScale)
    imp_nodes[[1]] <- impurete
  }

  id_feuille <- rep(1,length(Y_boot$id)) #### localisation des feuilles de l'arbre
  id_feuille_prime <- id_feuille
  feuilles_courantes <- unique(id_feuille)
  feuilles_terminales <- NULL

  for (p in 1:(length(unique(Y_boot$id))/2-1)){

    cat(paste0(p, "\n"))

    count_split <- 0
    for (i in 1:length(feuilles_courantes)){
      # Il faut que l'on regarde le tirage des variables de manière aléatoire :
      V <- NULL
      for (v in Inputs){
        V <- c(V, rep(get(v)$type,dim(get(v)$X)[length(dim(get(v)$X))]))
      }

      # mtry des espaces
      variables <- sample(V,mtry) # Maintenant on sait combien on doit en tirer dans chaque espace
      # On ne va regarder que les espaces tirés :
      split.spaces <- unique(variables)

      # variables <- sample(c(1:dim(X_boot$X[,,drop=FALSE])[2]),mtry)
      w <- which(id_feuille==feuilles_courantes[i])
      wXCurve <- NULL
      wXScalar <- NULL
      wXFactor <- NULL

      for (l in unique(Y_boot$id[w])){
        if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve_boot$id==l))
        if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar_boot$id==l))
        if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor_boot$id==l))
      }

      if (length(unique(Y_boot$id[w]))>1 & imp_nodes[[feuilles_courantes[i]]] >0){

        # On est ici

        # mtry des variables de chaque espace

        if (is.element("curve",split.spaces)==TRUE){
          tirageCurve <- sample(1:ncol(Curve$X),length(which(variables=="curve")))
          Curve_courant <- list(type = Curve_boot$type, X=Curve_boot$X[wXCurve,tirageCurve, drop=FALSE], id=Curve_boot$id[wXCurve, drop=FALSE], time=Curve_boot$time[wXCurve, drop=FALSE],
                                model = Curve_boot$model[tirageCurve])
        }

        if (is.element("scalar",split.spaces)==TRUE){
          tirageScalar <- sample(1:ncol(Scalar$X),length(which(variables=="scalar")))
          Scalar_courant <- list(type = Scalar_boot$type, X=Scalar_boot$X[wXScalar,tirageScalar, drop=FALSE], id=Scalar_boot$id[wXScalar, drop=FALSE])
        }

        if (is.element("factor",split.spaces)==TRUE){
          tirageFactor <- sample(1:ncol(Factor$X),length(which(variables=="factor")))
          Factor_courant <- list(type = Factor_boot$type, X=Factor_boot$X[wXFactor,tirageFactor, drop=FALSE], id=Factor_boot$id[wXFactor, drop=FALSE])
        }

        if (Y_boot$type=="curve"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w], id=Y_boot$id[w], time=Y_boot$time[w])
        }

        if (Y_boot$type=="surv"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w], id=Y_boot$id[w], comp=Y$comp)
        }

        if (Y_boot$type=="factor" || Y_boot$type=="scalar"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }


        F_SPLIT <- data.frame(TYPE = character(), Impurity = numeric(), stringsAsFactors = FALSE)
        decoupe <- 0

        # on test les meilleurs splits sur chacun des variables factor tire par mtry

        if (Y_courant$type == "surv"){
          #N_courant <- sum(Y_courant$Y[,2]!=0) # nb event
          N_courant <- sum(Y_courant$Y[,2]==cause)
        }else{
          N_courant <- length(unique(Y_courant$id))  # nb id
        }

        if (N_courant >= nodesize*2){ # si nb sujet/event inferieur a nodesize*2, feuille avec moins de nodesize

          if (is.element("factor",split.spaces)==TRUE){

            feuille_split_Factor <- var_split_MM(X = Factor_courant, Y = Y_courant, timeScale = timeScale,
                                                 nodesize = nodesize, cause = cause)

            if (feuille_split_Factor$Pure==FALSE){
              F_SPLIT <- merge(F_SPLIT,
                               data.frame(TYPE = "Factor", Impurity = feuille_split_Factor$impurete,
                                          stringsAsFactors = FALSE),
                               all = T)
              decoupe <- decoupe +1
            }
          }

          # on test les meilleurs splits sur chacun des markers tire par mtry

          if (is.element("curve",split.spaces)==TRUE){

            feuille_split_Curve <- var_split_MM(X = Curve_courant, Y = Y_courant, timeScale = timeScale,
                                                nsplit_option = nsplit_option, nodesize = nodesize,
                                                cause = cause)

            if (feuille_split_Curve$Pure==FALSE){
              F_SPLIT <- merge(F_SPLIT,
                               data.frame(TYPE = "Curve", Impurity = feuille_split_Curve$impurete,
                                          stringsAsFactors = FALSE),
                               all = T)
              decoupe <- decoupe +1
            }
          }

          if (is.element("scalar",split.spaces)==TRUE){

            feuille_split_Scalar <- var_split_MM(X = Scalar_courant, Y = Y_courant, timeScale = timeScale,
                                                 nsplit_option = nsplit_option, nodesize = nodesize,
                                                 cause = cause)

            if (feuille_split_Scalar$Pure==FALSE){
              F_SPLIT <- merge(F_SPLIT,
                               data.frame(TYPE = "Scalar", Impurity = feuille_split_Scalar$impurete,
                                          stringsAsFactors = FALSE),
                               all = T)
              decoupe <- decoupe +1
            }


          }

        }else{
          feuilles_terminales <- c(feuilles_terminales, feuilles_courantes[i])

          if (Y_courant$type=="surv"){
            Nevent <- sum(Y_courant$Y[,2]==cause) # nb event
          }else{
            Nevent <- NA
          }

          # add leafs to V_split
          V_split_node <- data.frame(type = "Leaf", num_noeud = feuilles_courantes[i], var_split = NA,
                                     var_summary = NA, threshold = NA, N = length(Y_courant$id),
                                     Nevent = Nevent, stringsAsFactors = FALSE)

          V_split <- merge(V_split, V_split_node, all = T)

          next()
        }

        if (decoupe>0){

          TYPE <- F_SPLIT[which.min(F_SPLIT[,2]),1]
          X <- get(TYPE)
          X_boot <- get(paste(TYPE,"_boot",sep=""))

          # on retrouve la repartition des individus OOB sur la variable qui minimise l'impurete

          feuille_split <- get(paste("feuille_split_",TYPE, sep=""))

          # on recupere la variable sur laquelle on a split

          vsplit_space <- get(paste("tirage",TYPE, sep=""))[feuille_split$variable]

          #if (imp_apres_split<imp_avant_split){

          gauche_id <- unique(Y_courant$id)[which(feuille_split$split==1)]
          droit_id <- unique(Y_courant$id)[which(feuille_split$split==2)]

          if (Y$type=="surv"){

            LNevent <- sum(Y_courant$Y[,2][which(feuille_split$split==1)]==cause)
            RNevent <- sum(Y_courant$Y[,2][which(feuille_split$split==2)]==cause)

            if (LNevent>=nodesize & RNevent>=nodesize){
              imp_nodes[[2*feuilles_courantes[i]]] <- Inf
              imp_nodes[[2*feuilles_courantes[i]+1]] <- Inf
            }else{
              feuilles_terminales <- c(feuilles_terminales, feuilles_courantes[i])
              Nevent <- sum(Y_courant$Y[,2]==cause) # nb event

              # add leafs to V_split
              V_split_node <- data.frame(type = "Leaf", num_noeud = feuilles_courantes[i], var_split = NA,
                                         var_summary = NA, threshold = NA, N = length(Y_courant$id),
                                         Nevent = Nevent, stringsAsFactors = FALSE)

              V_split <- merge(V_split, V_split_node, all = T)

              next()
            }

          }else{
            imp_nodes[[2*feuilles_courantes[i]]] <- feuille_split$impur_list[[1]]
            imp_nodes[[2*feuilles_courantes[i]+1]] <- feuille_split$impur_list[[2]]
          }

          if (Y_courant$type=="surv"){
            Nevent <- sum(Y_courant$Y[,2]==cause) # nb event
          }else{
            Nevent <- NA
          }

          # add node split to V_split
          V_split_node <- data.frame(type = TYPE, num_noeud = feuilles_courantes[i],
                                     var_split = vsplit_space, var_summary = feuille_split$variable_summary,
                                     threshold = feuille_split$threshold, N = length(Y_courant$id),
                                     Nevent = Nevent, stringsAsFactors = FALSE)

          V_split <- merge(V_split, V_split_node, all = T)

          model_param[[feuilles_courantes[i]]] <- feuille_split$model_param

          wY_gauche <- NULL
          wY_droit <- NULL
          w_gauche <- NULL
          w_droit <- NULL


          for (k in 1:length(gauche_id)){
            w_gauche <- c(w_gauche, which(X_boot$id==gauche_id[k]))
            wY_gauche <- c(wY_gauche, which(Y_boot$id==gauche_id[k]))
          }

          for (k in 1:length(droit_id)){
            w_droit <- c(w_droit, which(X_boot$id==droit_id[k]))
            wY_droit <- c(wY_droit, which(Y_boot$id==droit_id[k]))
          }

          id_feuille_prime[wY_gauche] <- 2*(feuilles_courantes[i])
          id_feuille_prime[wY_droit] <- 2*(feuilles_courantes[i])+1

          if (X$type=="curve"){
            trajG <- as.data.frame(cbind(X_boot$id[w_gauche], X_boot$time[w_gauche], X_boot$X[w_gauche,vsplit_space]))
            trajD <- as.data.frame(cbind(X_boot$id[w_droit], X_boot$time[w_droit], X_boot$X[w_droit,vsplit_space]))
            meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
            meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          }

          if (X$type=="factor"){
            meanFg <- unique(X_boot$X[w_gauche, vsplit_space])
            meanFd <- unique(X_boot$X[w_droit,vsplit_space])
          }

          if (X$type=="scalar"){
            meanFg <- mean(X_boot$X[w_gauche,vsplit_space])
            meanFd <- mean(X_boot$X[w_droit,vsplit_space])
          }


          hist_nodes[[2*(feuilles_courantes[i])]] <- meanFg
          hist_nodes[[2*(feuilles_courantes[i])+1]] <- meanFd
          count_split <- count_split+1

        }
      }
    }

    id_feuille <- id_feuille_prime
    feuilles_courantes <- setdiff(unique(id_feuille_prime), feuilles_terminales)

    if (count_split ==0 ){

      V_split <- V_split[order(V_split$num_noeud),]

      for (q in unique(id_feuille)){
        w <- which(id_feuille == q)
        if (Y$type=="curve"){
          Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
        }
        if(Y$type=="scalar"){
          Y_pred[[q]]<- mean(Y_boot$Y[w])
        }
        if(Y$type=="factor"){
          Table <- which.max(table(Y_boot$Y[w]))
          Y_pred[[q]] <-  as.factor(attributes(Table)$names)
        }
        if (Y$type=="surv"){
          datasurv <- data.frame(time_event = Y_boot$Y[w][,1], event = Y_boot$Y[w][,2])
          fit <- prodlim(Hist(time_event, event)~1, data = datasurv,
                         type = "risk")

          if (is.null(fit$cuminc)){
            pred <- data.frame(times=fit$time, traj=1-fit$surv) # 1-KM
          }else{
            pred <- data.frame(times=fit$time, traj=fit$cuminc[[cause]]) # CIF Aalen-Johansen
          }

          Y_pred[[q]] <- combine_times(pred = pred, newtimes = unique(Y$Y[,1]), type = "risk")
        }

      }
      if (Y$type=="factor"){
        Ylevels <- unique(Y_boot$Y)
        return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, boot=boot, Ylevels=Ylevels,
                    model_param = model_param))
      }
      return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, boot=boot,
                  model_param = model_param))
    }
  }

  V_split <- V_split[order(V_split$num_noeud),]

  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    if (Y$type=="curve"){
      Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
    }

    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y_boot$Y[w])
    }

    if(Y$type=="factor"){
      Table <- which.max(table(Y_boot$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }

    if (Y$type=="surv"){
      datasurv <- data.frame(time_event = Y_boot$Y[w][,1], event = Y_boot$Y[w][,2])
      fit <- prodlim(Hist(time_event, event)~1, data = datasurv,
                     type = "risk")

      if (is.null(fit$cuminc)){
        pred <- data.frame(times=fit$time, traj=1-fit$surv) # 1-KM
      }else{
        pred <- data.frame(times=fit$time, traj=fit$cuminc[[cause]]) # CIF Aalen-Johansen
      }

      Y_pred[[q]] <- combine_times(pred = pred, newtimes = unique(Y$Y[,1]), type = "risk")
    }

  }
  if (Y$type=="factor"){
    Ylevels <- unique(Y_boot$Y)
    return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, Ylevels=Ylevels, boot=boot))
  }
  return(list(feuilles = id_feuille,Ytype=Y_boot$type, idY=Y_boot$id, V_split=V_split, hist_nodes=hist_nodes, Y_pred= Y_pred, time=time, Y=Y, boot=boot))
}
