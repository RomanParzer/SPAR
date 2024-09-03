screening_ridge_lambda0 <- function(z, yz, family) {
  n <- NROW(z)
  p <- NCOL(z)
  tmp_sc <- apply(z,2,function(col)sqrt(var(col)*(n-1)/n))
  z2 <- scale(z,center=colMeans(z),scale=tmp_sc)
  lam_max <- 1000 * max(abs(t(yz)%*%z2[,tmp_sc>0]))/n*family$mu.eta(family$linkfun(mean(yz)))/family$variance(mean(yz))

  glmnet_res <- glmnet::glmnet(x=z, y=yz, family = family, alpha=0,lambda.min.ratio = min(0.01,1e-7/n / lam_max))
  lam <- min(glmnet_res$lambda)
  scr_coef <- coef(glmnet_res,s=lam)[-1]

  scr_coef
}

screening_marglik <- function(z, yz, family) {
  n <- nrow(z)
  scr_coef <- apply(z,2,function(zj){
    glm_res <- glm(yz~zj,family=family,start=c(1,0))
    glm_res$coefficients[2]
  })
  scr_coef
}

screening_corr <- function(z, yz, family) {
  n <- nrow(z)
  if (family$family != "gaussian") warning("Screening based on correlation is typically employed for Gaussian variables.")
  scr_coef <- sapply(seq_len(NCOL(z)), function(i) cor(yz, z[,i]))
  scr_coef
}

