#' Function to compute individual random effects using hlme output parameters
#'
#' @param model output object from hlme function
#' @param formula list of formula for fixed and random part
#' @param data data to compute random effect
#'
#' @importFrom splines ns
#'
#' @return
#' @export
#'
#' @import lcmm
#'
#' @examples
predRE <- function(model, formula, data){

  subject <- "id"
  dataNA <- na.omit(data)
  beta <- model$beta

  # Variance-covariance matrix of the random-effects

  B <- matrix(0, ncol = sum(model$idea0), nrow = sum(model$idea0))
  B[upper.tri(B,diag=TRUE)] <- model$varcov
  B <- t(B)
  B[upper.tri(B,diag=TRUE)] <- model$varcov

  se <- model$stderr^2 # residual variance error
  Z <- model.matrix(formula$random, dataNA) # random design matrix

  X <- model.matrix(formula$fixed, dataNA) # fixed design matrix

  Y <- model.matrix(reformulate(as.character(formula$fixed)[2], intercept = FALSE),
                    dataNA) # outcome

  bi <- matrix(NA, nrow = length(unique(data$id)), ncol = ncol(B),
               dimnames = list(unique(data$id), colnames(Z)))

  for (id in unique(data$id)){

    row.id <- which(dataNA$id==id)
    Zi <- Z[row.id, , drop = FALSE]
    Xi <- X[row.id, ]
    Yi <- Y[row.id, ]
    Vi <- Zi%*%B%*%t(Zi) + se*diag(length(row.id))
    b <- tryCatch(B%*%t(Zi)%*%solve(Vi)%*%(Yi-Xi%*%beta),
                  error = function(e) return(rep(NA, ncol(B)))) # solve issue

    bi[rownames(bi)==id,] <- b

  }

  return(list(bi = bi))

}
