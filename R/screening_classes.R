#' Constructor function for building screeningcoef objects
#'
#' Creates an object class screeningcoef using arguments passed by user.
#' @param name character
#' @param generate_scrcoef function for generating the screening coefficient. This
#' function should have with arguments \code{scrcoef}, which is a screeningcoef
#' object, and \code{data}, which is a list of x (matrix of predictors used as input in
#' \link{\code{spar}}) and y (vector of responses used in  \link{\code{spar}}).
#' @return a function which in turn creates a function which in turn creates an
# ' object of class \code{"screeningcoef"}
#' @description
#' The created function will return a object of class \code{"screeningcoef"} which
#' constitutes of a list. The attributes of the generating object will include by
#' default \code{type}, which can take one of two values \code{"prob"} (indicating
#' probabilistic screening should be employed) or
#' \code{"fixed"} (indicating that the top \code{nscreen} variables should be employed).
#' @export
constructor_screeningcoef <- function(name, generate_scrcoef) {
  ## Checks
  stopifnot(names(formals(generate_scrcoef)) %in% c("scrcoef", "data"))
  ## Function to return
  function(..., control = list()) {
    out <- list(name = name,
                generate_scrcoef = generate_scrcoef,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    if (is.null(attr(out, "split_data"))) attr(out, "split_data") <- FALSE
    if (is.null(attr(out, "type"))) {
      attr(out, "type") <- "prob"
    } else {
      stopifnot("'type' must be either 'prob' or 'fixed'." =
                  attr(out, "type") == "prob" | attr(out, "type") == "fixed")
    }
    class(out) <- c("screeningcoef")
    return(out)
  }
}
#' Function which generates the screening coefficients
#' @param scrcoef an object of class screeningcoef
#' @param data list of  \code{x} which is a matrix of covariates and
#'    \code{y} which is the response vector
#' @keywords internal
get_scrcoef <- function(scrcoef, data) {
  coef <- scrcoef$generate_scrcoef(scrcoef, data = data)
  coef
}

#'
#' Screening coefficient based on marginal GLMs
#'
generate_scrcoef_marglik <- function(scrcoef, data) {
  y <- data$y
  x <- data$x
  if (is.null(scrcoef$control$family)) {
    scrcoef$control$family <- attr(scrcoef, "family")
  }
  coefs <- apply(x, 2, function(xj){
    glm_res <- do.call(function(...) glm(y ~ xj,  ...),
                       scrcoef$control)
    glm_res$coefficients[2]
  })
  coefs
}
#' Screening coefficient based on marginal GLMs
#'
#' Creates an object class screeningcoef using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{screeningcoef}" object
#' @param control list of controls to be passed to the screening function
#' @return object of class "\code{screeningcoef}" which is a list with elements
#' \code{name} (character), \code{generate_scrcoef} (function), and \code{control}
#' (list of controls passes as an argument)
#' @description
#' No arguments need to be passed.
#'
#' @export
#'
screen_marglik <- constructor_screeningcoef(
  "screen_marglik",
  generate_scrcoef = generate_scrcoef_marglik)


#'
#' Screening coefficient based on correlation
#'
generate_scrcoef_corr <- function(scrcoef, data) {
  y <- data$y
  x <- data$x
  coefs <- apply(x, 2, function(xj) {
    do.call(function(...) cor(y, xj, ...),
            scrcoef$control)
  })
  coefs
}
#' Screening coefficient based  on correlation
#'
#' Creates an object class "\code{screeningcoef}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{screeningcoef}" object
#' @param control list of controls to be passed to the screening function
#' @return object of class "\code{screeningcoef}" which is a list with elements
#' \code{name} (character), \code{generate_scrcoef} (function), and \code{control}
#' (list of controls passes as an argument)
#' @description
#' No arguments need to be passed.
#'
#' @export
#'
screen_corr <- constructor_screeningcoef(
  "screen_corr",
  generate_scrcoef = generate_scrcoef_corr)
#'
#' Screening coefficient based  on ridge coefficients
#'
generate_scrcoef_ridge <- function(scrcoef, data) {
  z <- data$x
  yz <- data$y
  n <- NROW(z)
  p <- NCOL(z)
  if (is.null(scrcoef$control$family)) {
    scrcoef$control$family <- attr(scrcoef, "family")
  }
  family <- scrcoef$control$family
  # Set alpha to zero unless otherwise specified
  if (is.null(scrcoef$control$alpha)) scrcoef$control$alpha <- 0
  # Set lambda.min.ration to close to zero unless otherwise specified
  if (is.null(scrcoef$control$lambda.min.ratio)) {
    tmp_sc <- apply(z, 2, function(col) sqrt(var(col)*(n-1)/n))
    z2 <- scale(z, center = colMeans(z), scale = tmp_sc)
    ytZ <- crossprod(yz, z2[,tmp_sc > 0])
    lam_max <- 1000 * max(abs(ytZ))/n * family$mu.eta(family$linkfun(mean(yz)))/
      family$variance(mean(yz))
    scrcoef$control$lambda.min.ratio <- min(0.01, 1e-4 / lam_max)
  }
  # Obtain Ridge coefs GLMNET
  glmnet_res <- do.call(function(...) glmnet(x = z, y = yz, ...),
                        scrcoef$control)

  if (family$family == "gaussian") {
    dev.ratio_cutoff <- 0.999
  } else {
    dev.ratio_cutoff <- 0.8
  }
  lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= dev.ratio_cutoff])
  scr_coef <- coef(glmnet_res, s = lam)[-1]
  scr_coef
}

#' Screening coefficient based  on ridge coefficients
#'
#' Creates an object class "\code{screeningcoef}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{screeningcoef}" object
#' @param control list of controls to be passed to the screening function
#' @return object of class "\code{screeningcoef}" which is a list with elements
#' \code{name} (character), \code{generate_scrcoef} (function), and \code{control}
#' (list of controls passes as an argument)
#' @description
#' No arguments need to be passed.
#'
#' @export
#'
screen_ridge <- constructor_screeningcoef(
  "screen_ridge",
  generate_scrcoef = generate_scrcoef_ridge)


#' print.screeningcoef
#'
#' Print method for screeningcoef object
#' @param x description
#' @param ... further arguments passed to or from other methods
#' @return text summary
#'
#' @export
print.screeningcoef <- function(x, ...) {
  cat(paste0("Name: ", x$name), "\n")
  cat("Main attributes:", "\n")
  cat("split_data:",  attr(x, "split_data"), "\n")
  cat("type:",  attr(x, "type"), "\n")
  cat("number of screened variables:",  attr(x, "nscreen"), "\n")
  imp_vals <- attr(x, "importance")
  if (!is.null(imp_vals)) {
    cat("importance:",
        sprintf("num [1:%d] %s ...", length(imp_vals),
                paste(round(imp_vals[1:5], 3), collapse = " ")), "\n")
  }
}
