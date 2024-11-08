#' Constructor function for building randomprojection objects
#'
#' Creates an object class randomprojection using arguments passed by user.
#' @param name character
#' @param generate_fun function for generating the random projection matrix. This
#' function should have with arguments \code{rp}, which is a randomprojection
#' object, \code{m}, the target dimension and a vector of indexes
#' \code{included_vector} which shows the column index of the original variables in the
#' \code{x} matrix to be projected using the random projection. This is needed
#' due to the fact that screening can be employed pre-projection.
#' @param update_data_fun function for updating the randomprojection object with
#' information from the data. This
#' function should have with arguments \code{rp}, which is a randomprojection
#' object and \code{data}, which is a list containing `x` (the matrix of predictors)
#' and `y` the vector of responses.
#' @param update_rpm_w_data function for updating the random projection matrix with data.
#' This can be used for the case where a list of random projection matrices is
#' provided by argument \code{RPMs}. In this case, the random structure is kept
#' fixed, but the data-dependent part gets updated with the provided data. Defaults
#' to NULL. If not provided, the values of the provided RPMs do not change.
#' @param control list of controls for random projection. Can include minimum and
#' maximum dimension for the projection defaults to
#' \code{list(mslow = NULL, msup = NULL)}
#' @return a function which in turn creates an object of class randomprojection
#' @description
#' No arguments need to be passed. The entries of the matrix are generated from
#' a standard normal distribution.
#'
#' @export
constructor_rp <- function(name, generate_fun, update_data_fun = NULL,
                           update_rpm_w_data = NULL,
                           control = list(mslow = NULL, msup = NULL)) {
  ## Checks
  stopifnot(names(formals(generate_fun)) %in% c("rp", "m", "included_vector"))
  if (!is.null(update_data_fun)) {
    stopifnot(names(formals(update_data_fun)) %in% c("rp", "data"))
  }
  ## Function to return
  function(...) {
    out <- list(name = name,
                generate_rp_fun = generate_fun,
                update_data_fun = update_data_fun,
                control = control)
    attr <- list2(...)
    attributes(out) <- c(attributes(out), attr)
    if (is.null(attr(out, "data"))) attr(out, "data") <- FALSE
    class(out) <- c("randomprojection")
    return(out)
  }
}

#' Function which works on all random projection objects
#' @param rp an object of class randomprojection
#' @param m integer goal dimension of the projection
#' @param included_vector a vector containing column index of the original variables in the
#' \code{x} matrix to be projected using the random projection
#'
#' @keywords internal
get_rp <- function(rp, m, included_vector) {
  RM <- rp$generate_rp_fun(rp, m, included_vector)
  return(RM)
}

#' Gaussian random projection matrix
#'
#' Creates an object class randomprojection using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @return object of class randomprojection with is a list with elements name,
#' generate_rp_fun, update_data_fun, control
#' @description
#' No arguments need to be passed. The entries of the matrix are generated from
#' a standard normal distribution.
#'
#' @export
rp_gaussian <- function(...) {
  out <- list(name = "rp_gaussian",
              generate_rp_fun = generate_gaussian,
              control = list(mslow = NULL, msup = NULL))
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  attr(out, "data") <- FALSE
  class(out) <- c("randomprojection")
  return(out)
}
generate_gaussian <- function(rp, m, included_vector) {
  p <- length(included_vector)
  vals <- rnorm(m * p)
  RM <- matrix(vals, nrow = m, ncol = p)
  RM <- Matrix(RM, sparse = TRUE)
  return(RM)
}

#' Sparse random projection matrix
#'
#' Creates an object class randomprojection using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix. The possible argument is \code{psi} in (0,1] which determines
#' the level of sparsity in the matrix.
#' @return object of class randomprojection
#' @description
#' The sparse matrix used in \insertCite{ACHLIOPTAS2003JL}{SPAR} with entries equal to
#' \eqn{\Psi_{ij} = \pm 1/\sqrt{\psi}} with probability \eqn{\psi/2} and zero otherwise
#' for \eqn{\psi\in (0,1]}. Default is \code{psi = 1}.
#' @references{
#'   \insertRef{ACHLIOPTAS2003JL}{SPAR}
#' }
#' @export
rp_sparse <- function(...) {
  out <- list(name = "rp_sparse",
              generate_rp_fun = generate_sparse,
              control = list(mslow = NULL, msup = NULL))
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  if (is.null(attr(out, "psi"))) attr(out, "psi") <- 1
  attr(out, "data") <- FALSE
  class(out) <- c("randomprojection")
  return(out)
}

generate_sparse <- function(rp, m, included_vector) {
  p <- length(included_vector)
  psi <- attr(rp, "psi")
  if (psi > 1 | psi <= 0) stop("For a sparse rpm, psi should lie in interval (0,1].")
  v <- sample(c(-1, 0, 1), size = m * p,
              prob = c(psi/2, 1 - psi, psi/2), replace=TRUE)
  RM <- matrix(v/sqrt(psi), nrow = m, ncol = p)
  RM <- RM[rowSums(abs(RM)) > 0, ]
  RM <- Matrix(RM, sparse = TRUE)
  return(RM)
}

#' Sparse embedding matrix
#'
#' Creates an object class randomprojection using arguments passed by user.
#' @param ... includes arguments which can be passed as attributes to the random
#' projection matrix
#' @return object of class randomprojection
#' @description
#' The entries of the matrix are generated based on \insertCite{Clarkson2013LowRankApprox}{SPAR}.
#' @references{
#'   \insertRef{Clarkson2013LowRankApprox}{SPAR}
#' }
#' @export
rp_cw <- function(...) {
  out <- list(name = "rp_cw",
              generate_rp_fun = generate_cw,
              update_data_rp = update_data_cw,
              update_rpm_w_data = update_rpm_w_data_cw,
              control = list(mslow = NULL, msup = NULL))
  attr <- list2(...)
  attributes(out) <- c(attributes(out), attr)
  if (is.null(attr(out, "data"))) attr(out, "data") <- FALSE
  class(out) <- c("randomprojection")
  return(out)
}

generate_cw <- function(rp, m, included_vector) {
  p <- length(included_vector)
  use_data <- attr(rp, "data")
  if (!use_data) {
    diagvals <- sample(c(-1, 1), p, replace = TRUE)
  } else {
    if (is.null(attr(rp, "diagvals"))) stop("Must provide vector of coefficients for data-driven RP.")
    diagvals <- attr(rp, "diagvals")[included_vector]
  }
  goal_dims <- sample(m, p, replace = TRUE)
  counter <- 0
  # remove zero rows
  for (goal_dim in seq_len(m)) {
    if (sum(goal_dims==(goal_dim-counter))==0) {
      goal_dims[goal_dims > goal_dim - counter] <- goal_dims[goal_dims>goal_dim-counter]-1
      counter <- counter + 1
    }
  }
  RM <- Matrix(0, nrow = m - counter, ncol = p,sparse = TRUE)
  RM@i <- as.integer(goal_dims - 1)
  RM@p <- 0:p
  RM@x <- diagvals
  return(RM)
}

update_data_cw <- function(rp, data) {
  ## data should be list(x = x, y = y)
  if (is.null(data$x) || is.null(data$y)) stop("data must be a list of x, y.")
  z <- data$x
  yz <- data$y
  n <- NROW(z)
  p <- NCOL(z)
  family <- attr(rp, "family")
  tmp_sc <- apply(z,2,function(col)sqrt(var(col)*(n-1)/n))
  z2 <- scale(z,center=colMeans(z),scale=tmp_sc)
  ytZ <- crossprod(yz, z2[,tmp_sc>0])
  lam_max <- 1000 * max(abs(ytZ))/n*family$mu.eta(family$linkfun(mean(yz)))/family$variance(mean(yz))
  if (family$family=="gaussian") {
    dev.ratio_cutoff <- 0.999
  } else {
    dev.ratio_cutoff <- 0.8
  }
  glmnet_res <- glmnet(x=z, y=yz, family = family, alpha=0,
                       lambda.min.ratio = min(0.01,1e-4 / lam_max))
  lam <- min(glmnet_res$lambda[glmnet_res$dev.ratio<=dev.ratio_cutoff])
  scr_coef <- coef(glmnet_res,s=lam)[-1]
  inc_probs <- abs(scr_coef)
  max_inc_probs <- max(inc_probs)
  attr(rp, "diagvals") <- scr_coef/max_inc_probs
  return(rp)
}

update_rpm_w_data_cw <- function(rpm, rp, included_vector) {
  rpm@x <-  attr(rp, "diagvals")[included_vector]
  return(rpm)
}


