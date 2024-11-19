#' Constructor function for building screencoef objects
#'
#' Creates an object class screencoef using arguments passed by user.
#' @param name character
#' @param generate_fun function for generating the screening coefficient. This
#' function should  with arguments \code{object}, which is a "\code{screencoef}"
#' object, and \code{data}, which is a list of x (matrix of predictors used as input in
#' \link{\code{spar}}) and y (vector of responses used in  \link{\code{spar}}).
#' @return a function which in turn creates a function which in turn creates an
# ' object of class \code{"screencoef"}
#' @description
#' The created function will return a object of class "\code{screencoef}" which
#' constitutes of a list. The attributes of the generating object will include by
#' default \code{type}, which can take one of two values \code{"prob"} (indicating
#' probabilistic screening should be employed) or
#' \code{"fixed"} (indicating that the top \code{nscreen} variables should be employed).
#' @export
constructor_screencoef <- function(name, generate_fun) {
  ## Checks
  args_generate_fun <- formals(generate_fun)
  stopifnot("Function generate_fun should contain two arguments: an object
            of class \"screencoef\" and 'data'." =
              length(args_generate_fun) == 2)
  stopifnot("Function generate_fun should contain argument 'data'." =
              "data" %in% names(args_generate_fun))
  ## Function to return
  function(..., control = list()) {
    out <- list(name = name,
                generate_fun = generate_fun,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    if (is.null(attr(out, "type"))) {
      attr(out, "type") <- "prob"
    } else {
      stopifnot(
        "'type' must be either 'prob' or 'fixed'." =
          (attr(out, "type") == "prob" | attr(out, "type") == "fixed")
      )
    }
    class(out) <- c("screencoef")
    return(out)
  }
}
#' Function which generates the screening coefficients
#' @param object an object of class screencoef
#' @param data list of  \code{x} which is a matrix of covariates and
#'    \code{y} which is the response vector
#' @keywords internal
get_screencoef <- function(object, data) {
  coef <- object$generate_fun(object, data = data)
  coef
}

#'
#' Generate screening coefficient based  on marginal likelihood in univariate GLMs
#' @param object  "\code{screencoef}" object
#' @param data list of x and y
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_marglik <- function(object, data) {
  y <- data$y
  x <- data$x
  if (is.null(object$control$family)) {
    object$control$family <- attr(object, "family")
  }
  coefs <- apply(x, 2, function(xj){
    glm_res <- do.call(function(...) glm(y ~ xj,  ...),
                       object$control)
    glm_res$coefficients[2]
  })
  coefs
}
#' Screening coefficient based on marginal GLMs
#'
#' Creates an object class "\code{screencoef}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{screencoef}" object
#' @param control list of controls to be passed to the screening function
#' @return object of class "\code{screencoef}" which is a list with elements:
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  for generating the screening coefficient. This function should have arguments \code{object}, which is a "\code{screencoef}" object, and \code{data}, which is a list of two elements \code{x} and \code{y} containing the matrix of standardized predictors and the vector of (standardized for Gaussian) responses.
#' }
#'
#' @description
#' Relies on \link[stats]{glm}.
#'
#' @export
#'
screen_marglik <- constructor_screencoef(
  "screen_marglik",
  generate_fun = generate_scrcoef_marglik)


#'
#' Generate screening coefficient based  on correlation
#' @param object  "\code{screencoef}" object
#' @param data list of x and y
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_cor <- function(object, data) {
  y <- data$y
  x <- data$x
  coefs <- apply(x, 2, function(xj) {
    do.call(function(...) cor(y, xj, ...),
            object$control)
  })
  coefs
}
#' Screening coefficient based  on correlation
#'
#' Creates an object class "\code{screencoef}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{screencoef}" object
#' @param control list of controls to be passed to the screening function
#' @return object of class "\code{screencoef}" which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  for generating the screening coefficient. This function should have arguments \code{object}, which is a "\code{screencoef}" object, and \code{data}, which is a list of two elements \code{x} and \code{y} containing the matrix of standardized predictors and the vector of (standardized for Gaussian) responses.
#' }
#'
#' @description
#' Relies on \link[stats]{cor}.
#'
#' @export
#'
screen_cor <- constructor_screencoef(
  "screen_cor",
  generate_fun = generate_scrcoef_cor)

#'
#' Screening coefficient based  on glmnet coefficients
#' @param object  "\code{screencoef}" object
#' @param data list of x and y
#' @return vector of screening coefficients of length p
#' @keywords internal
generate_scrcoef_glmnet <- function(object, data) {
  z <- data$x
  yz <- data$y
  n <- NROW(z)
  p <- NCOL(z)
  control_glmnet <-
    object$control[names(object$control)  %in% names(formals(glmnet))]

  if (is.null(control_glmnet$family)) {
    control_glmnet$family <- attr(object, "family")
  }
  family <- control_glmnet$family

  # Set alpha to zero unless otherwise specified
  if (is.null(control_glmnet$alpha)) control_glmnet$alpha <- 0
  # Set lambda.min.ration to close to zero unless otherwise specified
  if (is.null(control_glmnet$lambda.min.ratio)) {
    tmp_sc <- apply(z, 2, function(col) sqrt(var(col)*(n-1)/n))
    z2 <- scale(z, center = colMeans(z), scale = tmp_sc)
    ytZ <- crossprod(yz, z2[,tmp_sc > 0])
    lam_max <- 1000 * max(abs(ytZ))/n * family$mu.eta(family$linkfun(mean(yz)))/
      family$variance(mean(yz))
    control_glmnet$lambda.min.ratio <- min(0.01, 1e-4 / lam_max)
  }
  # Obtain Ridge coefs GLMNET
  glmnet_res <- do.call(function(...) glmnet(x = z, y = yz, ...),
                        control_glmnet)

  if (family$family == "gaussian") {
    dev.ratio_cutoff <- 0.999
  } else {
    dev.ratio_cutoff <- 0.8
  }
  lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio <= dev.ratio_cutoff])
  scr_coef <- coef(glmnet_res, s = lam)[-1]
  scr_coef
}

#' Screening coefficient based  on glmnet coefficients
#'
#' Creates an object class "\code{screencoef}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{screencoef}" object
#' @param control list of controls to be passed to the screening function
#' @return object of class "\code{screencoef}" which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{generate_fun}  for generating the screening coefficient. This function should have arguments \code{object}, which is a "\code{screencoef}" object, and \code{data}, which is a list of two elements \code{x} and \code{y} containing the matrix of standardized predictors and the vector of (standardized for Gaussian) responses.
#' }
#'
#' @description
#' Relies on \link[glmnet]{glmnet}.
#'
#' @export
#'
screen_glmnet <- constructor_screencoef(
  "screen_glmnet",
  generate_fun = generate_scrcoef_glmnet)


#' print.screencoef
#'
#' Print method for a "\code{screencoef}" object
#' @param x description
#' @param ... further arguments passed to or from other methods
#' @return text summary
#'
#' @export
print.screencoef <- function(x, ...) {
  cat(paste0("Name: ", x$name), "\n")
  cat("Main attributes:", "\n")
  cat("* proportion of data used for screening:",
      ifelse(is.null(attr(x, "split_data_prop")),
             1, attr(x, "split_data_prop")), "\n")
  cat("* number of screened variables:",
      ifelse(is.null(attr(x, "nscreen")),
             "not provided, will default to 2n",
             attr(x, "nscreen")), "\n")
  cat("* type:",  ifelse(attr(x, "type") == "prob",
                         "probabilistic screening",
                         "screening top nscreen variables"), "\n")
  imp_vals <- attr(x, "importance")
  out_imp <-  ifelse(!is.null(imp_vals),
                     sprintf("num [1:%d] %s ...", length(imp_vals),
                             paste(round(imp_vals[1:5], 3),
                                   collapse = " ")),
                     "not yet computed from the data.")
  cat("* screening coefficients:", out_imp,  "\n")
}
