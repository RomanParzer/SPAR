#' Function which generates the screening coefficients
#' @param screeningcoef an object of class randomprojection
#' @param data list of  \code{x} which is a matrix of covariates and
#'    \code{y} which is the response vector
#' @keywords internal
get_scrcoef <- function(scrcoef, data) {
  coef <- scrcoef$generate_scrcoef(scrcoef, data = data)
  coef
}
#' Screening coefficient based on marginal GLMs
#'
#' Creates an object class screeningcoef using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @return object of class screeningcoef which is a list with elements name,
#' generate_scrcoef, control
#' @description
#' No arguments need to be passed.
#'
#' @export
screen_marglik <- function(...) {
  out <- list(name = "screen_marglik",
              generate_scrcoef = generate_scrcoef_marglik,
              control = list(nscreen = NULL, split_data = FALSE))
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
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
              generate_scrcoef = generate_scrcoef_corr,
              control = list(nscreen = NULL,
                             split_data = FALSE))
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("screeningcoef")
  return(out)
}

generate_scrcoef_corr <- function(scrcoef, data) {
  y <- data$y
  x <- data$x
  n <- nrow(x)
  method <- attr(scrcoef, "method")
  if (is.null(method)) method <- "pearson"
  coefs <- apply(x, 2, function(xj) cor(y, xj, method = method))
  coefs
}

#' @export
screen_ridge <- function(...) {
  out <- list(name = "screen_ridge",
              generate_scrcoef = generate_scrcoef_ridge,
              control = list(nscreen = NULL, split_data = FALSE))
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
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
  z2 <- scale(z, center = colMeans(z), scale = tmp_sc)
  ytZ <- crossprod(yz, z2[,tmp_sc > 0])
  lam_max <- 1000 * max(abs(ytZ))/n*family$mu.eta(family$linkfun(mean(yz)))/family$variance(mean(yz))
  if (family$family=="gaussian") {
    dev.ratio_cutoff <- 0.999
  } else {
    dev.ratio_cutoff <- 0.8
  }
  glmnet_res <- glmnet(x=z, y=yz, family = family,
                       alpha=0,lambda.min.ratio = min(0.01,1e-4 / lam_max))
  lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= dev.ratio_cutoff])
  scr_coef <- coef(glmnet_res, s = lam)[-1]
  scr_coef
}
