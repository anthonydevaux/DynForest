#' Randomized Dynamic tree
#'
#' @param Curve [list]:
#' @param Scalar [list]:
#' @param Factor [list]:
#' @param Shape [list]:
#' @param Image [list]:
#' @param Y [list]:
#' @param mtry [integer]:
#' @param aligned.shape [logical]:
#' @param timeScale [numeric]:
#' @param ... : option
#'
#' @import kmlShape
#' @import RiemBase
#' @import stringr
#' @import Evomorph
#' @import geomorph
#' @import survival
#' @importFrom splines ns
#'
#' @keywords internal
Rtmax <- function(Curve=NULL, Scalar=NULL, Factor=NULL, Shape=NULL, Image=NULL,Y,mtry,aligned.shape=FALSE, timeScale=0.1,
                  splitrule=NULL, nsplit_option=NULL, nodesize=1, ...){

  inputs <- read.Xarg(c(Curve,Scalar,Factor,Shape,Image))
  Inputs <- inputs

  for (k in 1:length(Inputs)){
    str_sub(Inputs[k],1,1) <- str_to_upper(str_sub(Inputs[k],1,1))
  }

  impurity_feuilles <- NULL
  V_split <- NULL
  hist_nodes <- list()
  model_param <- list()
  id_boot <- unique(sample(unique(Y$id), length(unique(Y$id)), replace=TRUE))
  boot <- id_boot
  decoupe <- 1

  wXCurve <- NULL
  wXScalar <- NULL
  wXFactor <- NULL
  wXShape <- NULL
  wXImage <- NULL
  wY <- NULL

  for (k in id_boot){
    wY <- c(wY, which(Y$id==k))
    if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve$id==k))
    if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar$id==k))
    if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor$id==k))
    if (is.element("shape",inputs)==TRUE) wXShape <- c(wXShape, which(Shape$id==k))
    if (is.element("image",inputs)==TRUE) wXImage <- c(wXImage, which(Image$id==k))
  }

  Y_pred <- list()

  if (is.element("curve",inputs)==TRUE) Curve_boot <- list(type=Curve$type,   X=Curve$X[wXCurve,, drop=FALSE], id= Curve$id[wXCurve], time = Curve$time[wXCurve],
                                                           model=Curve$model) ### bootstrap pour les courbes
  if (is.element("scalar",inputs)==TRUE) Scalar_boot <- list(type=Scalar$type,   X=Scalar$X[wXScalar,, drop=FALSE], id= Scalar$id[wXScalar]) ### bootstrap pour les courbes
  if (is.element("factor",inputs)==TRUE) Factor_boot <- list(type=Factor$type,   X=Factor$X[wXFactor,, drop=FALSE], id= Factor$id[wXFactor])
  if (is.element("shape",inputs)==TRUE) Shape_boot <- list(type=Shape$type,   X=Shape$X[,,wXShape, , drop=FALSE], id= Shape$id[wXShape])
  if (is.element("image",inputs)==TRUE) Image_boot <- list(type=Image$type,   X=Image$X[,,wXImage, , drop=FALSE], id= Image$id[wXImage])


  if (Y$type=="curve" || Y$type=="surv") {Y_boot <- list(type=Y$type,Y=Y$Y[wY], id=Y$id[wY], time=Y$time[wY])} ### idem pour Y
  if (Y$type=="image" || Y$type=="shape") {Y_boot <- list(type=Y$type, Y=Y$Y[,,wY], id=Y$id[wY])}
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

  for (p in 1:(length(unique(Y_boot$id))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){
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
      w <- which(id_feuille==unique(id_feuille)[i])
      wXCurve <- NULL
      wXScalar <- NULL
      wXFactor <- NULL
      wXShape <- NULL
      wXImage <- NULL

      for (l in unique(Y_boot$id[w])){
        if (is.element("curve",inputs)==TRUE) wXCurve <- c(wXCurve, which(Curve_boot$id==l))
        if (is.element("scalar",inputs)==TRUE) wXScalar <- c(wXScalar, which(Scalar_boot$id==l))
        if (is.element("factor",inputs)==TRUE) wXFactor <- c(wXFactor, which(Factor_boot$id==l))
        if (is.element("shape",inputs)==TRUE) wXShape <- c(wXShape, which(Shape_boot$id==l))
        if (is.element("image",inputs)==TRUE) wXImage <- c(wXImage, which(Image_boot$id==l))
      }

      if (length(unique(Y_boot$id[w]))>1 & imp_nodes[[unique(id_feuille)[i]]] >0){

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

        if (is.element("shape",split.spaces)==TRUE){
          tirageShape <- sample(1:dim(Shape$X)[length(dim(Shape$X))],length(which(variables=="shape")))
          Shape_courant <- list(type = Shape_boot$type, X=Shape_boot$X[,,wXShape,tirageShape, drop=FALSE], id=Shape_boot$id[wXShape, drop=FALSE])
        }

        if (is.element("image",split.spaces)==TRUE){
          tirageImage <- sample(1:dim(Image$X)[length(dim(Image$X))],length(which(variables=="image")))
          Image_courant <- list(type = Image_boot$type, X=Image_boot$X[,,wXImage,tirageImage, drop=FALSE], id=Image_boot$id[wXImage])
        }

        if (Y_boot$type=="curve" || Y_boot$type=="surv"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w], id=Y_boot$id[w], time=Y_boot$time[w])
        }

        if (Y_boot$type=="image" || Y_boot$type=="shape"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[,,w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }


        if (Y_boot$type=="factor" || Y_boot$type=="scalar"){
          Y_courant <- list(type=Y_boot$type, Y=Y_boot$Y[w, drop=FALSE], id=Y_boot$id[w, drop=FALSE])
        }


        F_SPLIT <- NULL
        decoupe <- 0

        # on test les meilleurs splits sur chacun des variables factor tire par mtry

        if (length(unique(Y_courant$id)) >= nodesize*2){ # si nb sujet inferieur a nodesize*2, feuille avec moins de nodesize

          if (is.element("factor",split.spaces)==TRUE){

              if (splitrule=="MM"){
                feuille_split_Factor <- var_split_MM(Factor_courant,Y_courant,timeScale,
                                                          nodesize)
              }

            if (feuille_split_Factor$Pure==FALSE){
              F_SPLIT <- rbind(F_SPLIT,c("Factor",feuille_split_Factor$impurete))
              decoupe <- decoupe +1
            }
          }

          # on test les meilleurs splits sur chacun des markers tire par mtry

          if (is.element("curve",split.spaces)==TRUE){

              if (splitrule=="MM"){
                feuille_split_Curve <- var_split_MM(Curve_courant,Y_courant,timeScale,
                                                         nsplit_option,nodesize)
              }

            if (feuille_split_Curve$Pure==FALSE){
              F_SPLIT <- rbind(F_SPLIT,c("Curve",feuille_split_Curve$impurete))
              decoupe <- decoupe +1
            }
          }

          if (is.element("scalar",split.spaces)==TRUE){

              if (splitrule=="MM"){
                feuille_split_Scalar <- var_split_MM(Scalar_courant,Y_courant,timeScale,
                                                          nsplit_option,nodesize)
              }

            if (feuille_split_Scalar$Pure==FALSE){
              F_SPLIT <- rbind(F_SPLIT,c("Scalar",feuille_split_Scalar$impurete))
              decoupe <- decoupe +1
            }


          }

          if (is.element("shape",split.spaces)==TRUE){

            feuille_split_Shape <- ERvar_split(Shape_courant,Y_courant,timeScale, ntry=ntry)


            if (feuille_split_Shape$Pure==FALSE){
              F_SPLIT <- rbind(F_SPLIT,c("Shape",feuille_split_Shape$impurete))
              decoupe <- decoupe +1
            }

          }


          if (is.element("image",split.spaces)==TRUE){

            feuille_split_Image <- ERvar_split(Image_courant,Y_courant,timeScale, ntry=ntry)

            if (feuille_split_Image$Pure==FALSE){
              F_SPLIT <- rbind(F_SPLIT,c("Image",feuille_split_Image$impurete))
              decoupe <- decoupe +1
            }

          }

        }else{

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

          gauche_id <- unique(Y_boot$id[w])[which(feuille_split$split==1)]
          droit_id <- unique(Y_boot$id[w])[which(feuille_split$split==2)]

          if (Y$type=="surv"){
            imp_nodes[[2*unique(id_feuille)[i]]] <- Inf
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- Inf
          }
          else {
            imp_nodes[[2*unique(id_feuille)[i]]] <- feuille_split$impur_list[[1]]
            imp_nodes[[2*unique(id_feuille)[i]+1]] <- feuille_split$impur_list[[2]]
          }

          # split sur quel espace, quel noeud, quelle variable
          #V_split <- rbind(V_split,c(TYPE,unique(id_feuille)[i],vsplit_space))

          # split sur quel espace, quel noeud, quelle variable, quel resume, quel threshold
          V_split <- rbind(V_split,c(TYPE,unique(id_feuille)[i],vsplit_space,
                                     feuille_split$variable_summary, feuille_split$threshold))

          model_param[[unique(id_feuille)[i]]] <- feuille_split$model_param

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


          id_feuille_prime[wY_gauche] <- 2*(unique(id_feuille)[i])
          id_feuille_prime[wY_droit] <- 2*(unique(id_feuille)[i])+1

          #print(paste("Split on the variable", vsplit_space, "on the space of ", paste(TYPE,"s",sep="")))

          if (X$type=="curve"){
            trajG <- as.data.frame(cbind(X_boot$id[w_gauche], X_boot$time[w_gauche], X_boot$X[w_gauche,vsplit_space]))
            trajD <- as.data.frame(cbind(X_boot$id[w_droit], X_boot$time[w_droit], X_boot$X[w_droit,vsplit_space]))
            meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
            meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          }

          if (X$type=="shape"){
            meanFg <- mshape(X_boot$X[,,w_gauche,vsplit_space])
            meanFd <- mshape(X_boot$X[,,w_droit,vsplit_space])
          }

          if (X$type=="image"){
            meanFg <- rbase.mean(riemfactory(X_boot$X[,,w_gauche,vsplit_space]))
            meanFd <- rbase.mean(riemfactory(X_boot$X[,,w_droit,vsplit_space]))
            hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg$x
            hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd$x
          }

          if (X$type=="factor"){
            meanFg <- unique(X_boot$X[w_gauche, vsplit_space])
            meanFd <- unique(X_boot$X[w_droit,vsplit_space])
          }

          if (X$type=="scalar"){
            meanFg <- mean(X_boot$X[w_gauche,vsplit_space])
            meanFd <- mean(X_boot$X[w_droit,vsplit_space])
          }



          hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
          hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
          count_split <- count_split+1

          feuilles_courantes <- unique(id_feuille_prime)
        }


      }
    }

    id_feuille <- id_feuille_prime

    if (count_split ==0 ){

      V_split <- data.frame(V_split)
      names(V_split) <- c("type","num_noeud", "var_split","var_summary","threshold")

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

        if (Y$type=="shape"){
          Y_pred[[q]] <-  mshape(Y_boot$Y[,,w, drop=FALSE])
        }

        if (Y$type=="image"){
          donnees <- riemfactory(Y_boot$Y[,,w, drop=FALSE])
          Y_pred[[q]] <- rbase.mean(donnees)$x
        }

        if (Y$type=="surv"){
          donnees = survfit(Surv(Y_boot$time[w], Y_boot$Y[w])~1)
          # CHF
          Y_pred[[q]] <- data.frame(times=donnees$time, traj=donnees$cumhaz)
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

  V_split <- data.frame(V_split)
  names(V_split) <- c("type","num_noeud", "var_split","var_summary","threshold")

  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    if (Y$type=="curve"){
      Y_pred[[q]] <- kmlShape::meanFrechet(data.frame(Y_boot$id[w], Y_boot$time[w], Y_boot$Y[w]))
    }

    if (Y$type=="image"){
      donnees <- riemfactory(Y_boot$Y[,,w, drop=FALSE])
      Y_pred[[q]] <- rbase.mean(donnees)$x
    }

    if(Y$type=="scalar"){
      Y_pred[[q]]<- mean(Y_boot$Y[w])
    }

    if(Y$type=="factor"){
      Table <- which.max(table(Y_boot$Y[w]))
      Y_pred[[q]] <-  as.factor(attributes(Table)$names)
    }

    if (Y$type=="shape"){
      Y_pred[[q]] <- mshape(Y_boot$Y[,,w, drop=FALSE])
    }

    if (Y$type=="surv"){
      donnees = survfit(Surv(Y_boot$time[w], Y_boot$Y[w])~1)
      # CHF
      Y_pred[[q]] <- data.frame(times=donnees$time, traj=donnees$cumhaz)
    }

  }
  if (Y$type=="factor"){
    Ylevels <- unique(Y_boot$Y)
    return(list(feuilles = id_feuille, idY=Y_boot$id,Ytype=Y_boot$type, V_split=V_split, hist_nodes=hist_nodes, Y_pred = Y_pred, time = time, Y=Y, Ylevels=Ylevels, boot=boot))
  }
  return(list(feuilles = id_feuille,Ytype=Y_boot$type, idY=Y_boot$id, V_split=V_split, hist_nodes=hist_nodes, Y_pred= Y_pred, time=time, Y=Y, boot=boot))
}
