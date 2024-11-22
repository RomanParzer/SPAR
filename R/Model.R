#' Model object for estimating penalized glms in marginal models in the ensemble
#'
#' Creates an object class "\code{sparmodel}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{sparmodel}" object
#' @param control list of controls to be passed to the model function
#' @return object of class "\code{sparmodel}" which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{model_fun}  for generating the screening coefficient.
#'   This function should have arguments \code{y}, vector of standardized responses,
#'   \code{z}, a matrix of projected predictors in each marginal model, and \code{object}, which is a "\code{sparmodel}" object. Returns a list with two elements: \code{gammas} which is the vector of regression coefficients for the projected predictors and \code{intercept} which is the intercept of the model.
#'  \item \code{update_sparmodel}  optional function for updating the sparmodel object. before the
#' start of the algorithm.
#' }
#' @description
#' Relies on \link[glmnet]{glmnet}.
#'
#' @export
#'
spar_glmnet <- function(..., control = list()) {
  out <-  list(name = "glmnet",
               model_fun = model_glmnet,
               update_sparmodel = update_sparmodel_glmnet,
               control = control)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("sparmodel")
  out
}

ols_fun <- function(y, z) {
  solve(crossprod(z), crossprod(z,y))
}
ols_fun_corrected <- function(y, z) {
  solve(crossprod(z) + 0.01*diag(ncol(z)), crossprod(z,y))
}

update_sparmodel_glmnet <- function(object) {
  family <- attr(object, "family")
  if (family$family=="gaussian" & family$link=="identity") {
    fit_family <- "gaussian"
  } else {
    if (family$family=="binomial" & family$link=="logit") {
      fit_family <- "binomial"
    } else if (family$family=="poisson" & family$link=="log") {
      fit_family <- "poisson"
    } else {
      fit_family <- family
    }
  }
  attr(object, "family") <- family
  object
}

model_glmnet <- function(y, z, object) {
  ## y - vector of n responses
  ## z - matrix with n rows
  if (is.null(object$control$family)) {
    object$control$family <- attr(object, "family")
  }
  if (is.null(object$control$alpha)) {
    object$control$alpha <- 0
  }
  family <- object$control$family

  if (family$family=="gaussian" & family$link=="identity") {
    gammas <- tryCatch(ols_fun(y, z),
                         error=function(error_message) {
                           return(ols_fun_corrected(y, z))
                         })
    intercept <- 0
  } else {
    glmnet_res <- do.call(function(...) glmnet(x = z, y = y, ...),
                          object$control)
    mar_coef <- coef(glmnet_res, s = min(glmnet_res$lambda))
    intercept <- mar_coef[1]
    gammas <- mar_coef[-1]
  }
  list(gammas = gammas, intercept = intercept)
}

#' Model object for estimating penalized glms in marginal models in the ensemble models in the ensemble
#'
#' Creates an object class "\code{sparmodel}" using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the
#' "\code{sparmodel}" object
#' @param control list of controls to be passed to the model function
#' @return object of class "\code{sparmodel}" which is a list with elements
#'
#' \itemize{
#'  \item \code{name} (character)
#'  \item \code{control} (list of controls passed as an argument)
#'  \item \code{model_fun}  for generating the screening coefficient.
#'   This function should have arguments \code{y}, vector of standardized responses,
#'   \code{z}, a matrix of projected predictors in each marginal model, and \code{object}, which is a "\code{sparmodel}" object. Returns a list with two elements: \code{gammas} which is the vector of regression coefficients for the projected predictors and \code{intercept} which is the intercept of the model
#' }
#' @description
#' Relies on \link[stats]{glm.fit}.
#'
#' @export
#'
spar_glm <- function(..., control = list()) {
  out <-  list(name = "glm",
               model_fun = model_glm,
               control = control)
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  class(out) <- c("sparmodel")
  out
}

model_glm <- function(y, z, object) {
  ## y - vector of n responses
  ## z - matrix with n rows
  if (is.null(object$control$family)) {
    object$control$family <- attr(object, "family")
  }
  glm_res <- do.call(function(...) glm.fit(x = z, y = y, ...),
                          object$control)
  intercept <- coef(res)[1]
  gammas <- coef(res)[-1]
  list(gammas = gammas, intercept = intercept)
}

