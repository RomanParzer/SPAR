get_scrcoef <- function(scrcoef, data) {
  coef <- scrcoef$generate_scrcoef(scrcoef, data = data)
  coef
}
#' @export
screen_marglik <- function(...) {
  out <- list(name = "screen_marglik",
              generate_scrcoef = generate_scrcoef_marglik)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  if (is.null(attr(out, "family"))) {
    warning("No family provided for screening coefficient based on marginal likelihood in glms. Using gaussian() by default.")
    attr(out, "family") <- gaussian()
  }
  class(out) <- c("screeningcoef")
  return(out)
}
generate_scrcoef_marglik <- function(scrcoef, data) {
  y <- data$y
  x <- data$x
  n <- nrow(x)
  family <- attr(scrcoef, "family")
  coefs <- apply(x, 2, function(xj){
    glm_res <- glm(y ~ xj, family = family,
                   start = c(1, 0))
    glm_res$coefficients[2]
  })
  coefs
}

#' @export
screen_corr <- function(...) {
  out <- list(name = "screen_corr",
              generate_scrcoef = generate_scrcoef_corr)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("screeningcoef")
  return(out)
}
generate_scrcoef_corr <- function(scrcoef, data) {
  y <- data$y
  x <- data$x
  n <- nrow(x)
  family <- attr(scrcoef, "family")
  coefs <- apply(x, 2, function(xj) cor(y, xj, ))
  coefs
}

#' @export
screen_ridge <- function(...) {
  out <- list(name = "screen_ridge",
              generate_scrcoef = generate_scrcoef_ridge)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  if (is.null(attr(out, "family"))) {
    warning("No family provided for screening coefficient based on marginal likelihood in glms. Using gaussian() by default.")
    attr(out, "family") <- gaussian()
  }
  class(out) <- c("screeningcoef")
  return(out)
}
generate_scrcoef_ridge <- function(scrcoef, data) {
  z <- data$x
  yz <- data$y
  family <- attr(scrcoef, "family")
  n <- NROW(z)
  p <- NCOL(z)
  tmp_sc <- apply(z, 2, function(col) sqrt(var(col)*(n-1)/n))
  z2 <- scale(z,center = colMeans(z),scale = tmp_sc)
  ytZ <- crossprod(yz, z2[,tmp_sc>0])
  lam_max <- 1000 * max(abs(ytZ))/n*family$mu.eta(family$linkfun(mean(yz)))/family$variance(mean(yz))
  if (family$family=="gaussian") {
    dev.ratio_cutoff <- 0.999
  } else {
    dev.ratio_cutoff <- 0.8
  }
  glmnet_res <- glmnet(x=z, y=yz, family = family,
                       alpha=0,lambda.min.ratio = min(0.01,1e-4 / lam_max))
  lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio<=dev.ratio_cutoff])
  scr_coef <- coef(glmnet_res,s=lam)[-1]
  scr_coef
}
