## 1D interpolation function ====================================================
Interpol1D <- function(y, t, tNew){
  # Function that interpolates 1d function from t grid to tNew grid

  yNew <- rep(NA, length(tNew))
  idx_tNew <- sapply(tNew, function(x) sum(x>=t))
  for (j in seq(length(tNew))){
    tj <- tNew[j]
    if (tj %in% t){
      yNew[j] <- y[idx_tNew[j]]
    } else {
      idx_tauj <- idx_tNew[j]
      tauj <- t[idx_tauj]
      taujplus <- t[idx_tauj+1]
      yNew[j] <- (taujplus - tj)/(taujplus - tauj)*y[idx_tauj] + (tj - tauj)/(taujplus-tauj)*y[idx_tauj+1]
    }
  }
  return(yNew)
}

## Covariance matrix interpolation function ====================================
InterpolCovMat <- function(CovMat, tMat, tNew){
  # Function that interpolates CovMat from tMat grid to tNew grid
  # add pre-treatment to remove tNew not in the range of tMat
  mat_tNew <- matrix(NA,nrow=length(tNew), ncol=length(tNew))
  idx_tNew <- sapply(tNew, function(x) sum(x>=tMat))

  for (j in seq(length(tNew))){

    tj <- tNew[j]
    idx_tauj <- idx_tNew[j]
    tauj <- tMat[idx_tauj]
    taujplus <- tMat[idx_tauj+1]
    # on parcourt la triangulaire inférieure
    for (k in seq(j)){
      tk <- tNew[k]
      idx_tauk <- idx_tNew[k]
      tauk <- tMat[idx_tauk]
      taukplus <-  tMat[idx_tauk+1]
      # matrix interpolation
      Gjk <- CovMat[idx_tauj,idx_tauk] +
        (tj-tauj)/(taujplus-tauj)*(CovMat[idx_tauj+1, idx_tauk]-CovMat[idx_tauj,idx_tauk]) +
        (tj-tauj)*(tk-tauk)/((taujplus-tauj)*(taukplus-tauk))*
        (CovMat[idx_tauj+1, idx_tauk+1] - CovMat[idx_tauj+1, idx_tauk] - CovMat[idx_tauj,idx_tauk+1] + CovMat[idx_tauj,idx_tauk]) +
        (tk-tauk)/(taukplus-tauk)*(CovMat[idx_tauj,idx_tauk+1]-CovMat[idx_tauj,idx_tauk])
      if (j != k) {
        mat_tNew[j,k] <- mat_tNew[k,j] <- Gjk
      } else {
        mat_tNew[j,j] <- Gjk
      }
    }
  }
  return(mat_tNew)
}

pred_fpca_manual <- function(FPCAobj, dt_Ly_test, dt_Lt_test, dt_Lt_train){

  scores <- matrix(NA, nrow = length(dt_Ly_test), ncol = FPCAobj$selectK)
  min_dt_Lt_train <- min(unlist(dt_Lt_train), na.rm = TRUE)
  max_dt_Lt_train <- max(unlist(dt_Lt_train), na.rm = TRUE)

  for (i in seq(length(dt_Ly_test))){

    Lti <- dt_Lt_test[[i]][dt_Lt_test[[i]]>=min_dt_Lt_train & dt_Lt_test[[i]]<max_dt_Lt_train]
    Lyi <- dt_Ly_test[[i]][dt_Lt_test[[i]]>=min_dt_Lt_train & dt_Lt_test[[i]]<max_dt_Lt_train]
    Lti <- Lti[!is.na(Lyi)]; Lyi <- Lyi[!is.na(Lyi)];

    if (length(Lti) > 0){

      interpol_mui <- Interpol1D(FPCAobj$mu, FPCAobj$workGrid, Lti)
      interpol_FPCi <- apply(FPCAobj$phi, 2, function(x) return(Interpol1D(x, FPCAobj$workGrid, Lti)))
      interpolCov <- InterpolCovMat(FPCAobj$fittedCov, FPCAobj$workGrid, Lti)

      scores[i,] <- t(matrix(rep(FPCAobj$lambda, length(Lti)), nrow = length(Lti), byrow = TRUE) * interpol_FPCi) %*%
        solve(interpolCov + FPCAobj$sigma2*diag(length(Lti))) %*% (Lyi - interpol_mui)
    }
  }
  return(scores)
}

pred_fpca_manual2 <- function(workgrid, K, mu, FPCs, Cov, sigma2, lambda, min_dt_Lt_train, max_dt_Lt_train, dt_Ly_test, dt_Lt_test){

  scores <- matrix(NA, nrow = length(dt_Ly_test), ncol = K)
  for (i in seq(length(dt_Ly_test))){

    Lti <- dt_Lt_test[[i]][dt_Lt_test[[i]]>=min_dt_Lt_train & dt_Lt_test[[i]]<=max_dt_Lt_train]
    Lyi <- dt_Ly_test[[i]][dt_Lt_test[[i]]>=min_dt_Lt_train & dt_Lt_test[[i]]<=max_dt_Lt_train]
    Lti <- Lti[!is.na(Lyi)]; Lyi <- Lyi[!is.na(Lyi)];

    interpol_mui <- Interpol1D(mu, workgrid, Lti)
    interpol_FPCi <- apply(FPCs, 2, function(x) return(Interpol1D(x, workgrid, Lti)))

    interpolCov <- InterpolCovMat(Cov, workgrid, Lti)
    scores[i,] <- t(matrix(rep(lambda, length(Lti)), nrow = length(Lti), byrow = TRUE) * interpol_FPCi) %*%
      solve(interpolCov + sigma2*diag(length(Lti))) %*% (Lyi - interpol_mui)

  }

  return(scores)
}

pred_fdapace <- function(FPCAobj, dt_Ly_test, dt_Lt_test, dt_Lt_train){

  # censoring
  min_dt_Lt_train <- min(unlist(dt_Lt_train), na.rm = TRUE)
  max_dt_Lt_train <- max(unlist(dt_Lt_train), na.rm = TRUE)
  idx_uncensored <- lapply(dt_Lt_test, function(x) return(which(x>=min_dt_Lt_train & x<=max_dt_Lt_train)))
  dt_Lt_test_cens <- Map(function(a,b) return(a[b]), dt_Lt_test, idx_uncensored)
  dt_Ly_test_cens <- Map(function(a,b) return(a[b]), dt_Ly_test, idx_uncensored)

  pred_fpca_fdapace_test <- predict(FPCAobj, dt_Ly_test_cens, dt_Lt_test_cens)

  return(pred_fpca_fdapace_test$scores)
}
