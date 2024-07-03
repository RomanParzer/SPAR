screening_ridge_lambda0 <- function(z, yz, family) {
  n <- NROW(z)
  p <- NCOL(z)
  if (family == "gaussian") {
    if (p < n/2) {
      scr_coef <- tryCatch(
        solve(crossprod(z),crossprod(z, yz)),
        error=function(error_message) {
          return(solve(crossprod(z)+(sqrt(p)+sqrt(n)) * diag(p),crossprod(z,yz)))
        })
    } else if (p < 2*n) {
      scr_coef <- crossprod(z,solve(tcrossprod(z)+(sqrt(p)+sqrt(n)) * diag(n), yz))
    } else {
      solve_res <- tryCatch( solve(tcrossprod(z),yz),
                             error=function(error_message) {
                               return(solve(tcrossprod(z)+(sqrt(p)+sqrt(n))*diag(n),yz))
                             })
      scr_coef <- crossprod(z,solve_res)
    }
  } else {
    glmnet_res <- glmnet::glmnet(x=z, y=yz, family = family, alpha=0)
    lam <- min(glmnet_res$lambda)
    scr_coef <- coef(glmnet_res,s=lam)[-1]
  }
  scr_coef
}

screening_marglik <- function(z, yz, family) {
  n <- nrow(z)
  scr_coef <- sapply(seq_len(NCOL(z)), function(i)
    logLik(glm(formula = yz ~ z[,i], family = family)))
  scr_coef
}

screening_corr <- function(z, yz, family) {
  n <- nrow(z)
  if (family$family != "gaussian") warning("Screening based on correlation is typically employed for Gaussian variables.")
  scr_coef <- sapply(seq_len(NCOL(z)), function(i) cor(yz, z[,i]))
  scr_coef
}

