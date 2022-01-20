#' Internal function for simulations
#'
#' @importFrom lme4 ranef fixef VarCorr
#'
#' @export
#'
#' @examples
#'
predRE_simu <- function(model, data){

  subject <- "id"

  row_subject <- rownames(data)

  beta <- model$params$beta # fixed effect

  len.sd.re <- length(model$params$sd.re)

  if (len.sd.re==1){

    B <- as.matrix(model$params$sd.re)

  }else if (len.sd.re==2){

    covar <- model$params$sd.re[1] * model$params$sd.re[2] * model$params$cor.re

    B <- matrix(c(model$params$sd.re[1]^2, covar,
                  covar, model$params$sd.re[2]^2), ncol = 2)

  }else if (len.sd.re==3){

    covar01 <- model$params$sd.re[1] * model$params$sd.re[2] * model$params$cor.re[1]
    covar02 <- model$params$sd.re[1] * model$params$sd.re[3] * model$params$cor.re[2]
    covar12 <- model$params$sd.re[2] * model$params$sd.re[3] * model$params$cor.re[3]

    B <- matrix(c(model$params$sd.re[1]^2, covar01, covar02,
                  covar01, model$params$sd.re[2]^2, covar12,
                  covar02, covar12, model$params$sd.re[3]^2), ncol = 3)

  }else if (len.sd.re==4){

    B <- matrix(c(model$params$sd.re[1]^2, 0, 0, 0,
                  0, model$params$sd.re[2]^2, 0, 0,
                  0, 0, model$params$sd.re[3]^2, 0,
                  0, 0, 0, model$params$sd.re[4]^2), ncol = 4)

  }

  se <- model$params$sigmae^2 # residual variance error

  # random design matrix

  Z.formula <- model$random
  Z <- model.matrix(Z.formula, data)

  row_subject <- intersect(row_subject, rownames(Z))

  # fixed design matrix

  X.formula <- model$fixed
  X <- model.matrix(X.formula, data)

  row_subject <- intersect(row_subject, rownames(X))

  # outcome
  Y <- na.omit(data[,as.character(model$fixed)[2], drop = FALSE])

  row_subject <- intersect(row_subject, rownames(Y))

  data <- data[row_subject,]

  predRE <- matrix(NA, nrow = length(unique(data[,subject])), ncol = ncol(B),
                   dimnames = list(unique(data[,subject]), colnames(Z)))
  predRE_row <- 1

  for (ind_subject in unique(data[,subject])){

    ind <- rownames(data[which(data[, subject] == ind_subject),])

    Lsubject <- nrow(data[ind,])

    Z_i <- matrix(Z[ind,], nrow = Lsubject)
    V_i <- Z_i%*%B%*%t(Z_i) + se*diag(Lsubject)
    Y_i <- Y[ind,]
    X_i <- X[ind,]
    b_i <- B%*%t(Z_i)%*%solve(V_i)%*%(Y_i-X_i%*%beta)

    predRE[predRE_row,] <- b_i

    predRE_row <- predRE_row + 1
  }

  return(list(b_i = predRE,
              beta = beta,
              group = subject))

}
