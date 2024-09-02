
#' Generate CW Random Projection Matrix (Definition 2.2. in Parzer, Vana_Gür, Filzmoser 2023)
#'
#' @param m goal dimension
#' @param p initial dimension
#' @param coef p-vector of diagonal elements, random +/-1 by default
#' @returns sparse mxp matrix of class "dgCMatrix" (Matrix package)
#' @keywords internal
generate_cw_rp <- function(type.rpm, m, p, coef = NULL) {

  goal_dims <- sample(m, p, replace = TRUE)
  counter <- 0
  # remove zero rows
  for (goal_dim in seq_len(m)) {
    if (sum(goal_dims==(goal_dim-counter))==0) {
      goal_dims[goal_dims>goal_dim-counter] <- goal_dims[goal_dims>goal_dim-counter]-1
      counter <- counter + 1
    }
  }

  RM <- Matrix::Matrix(0, nrow = m - counter, ncol = p,sparse = TRUE)
  RM@i <- as.integer(goal_dims - 1)
  RM@p <- 0:p
  RM@x <- coef

  return(RM)
}

#' Sparse Random Projection Matrix
#'
#' @param m goal dimension
#' @param p dimension of predictors to be projected
#' @param psi a vector of length one representing the parameter of the random projection 0 < vals <= 1
#'
#' @return an m x p matrix of class "dgCMatrix"
#' @export
#'
#' @examples
#' generate_sparse_rp(10, 20)
generate_sparse_rp <- function(m, p, psi = NULL) {
  if (is.null(psi)) psi <- 0.1
  if (psi > 1 | psi <= 0) stop("For a sparse rpm, psi should lie in interval (0,1].")
  v <- sample(c(-1, 0, 1), size = m * p,
              prob = c(psi/2, 1 - psi, psi/2), replace=TRUE)
  RM <- matrix(v/sqrt(psi), nrow = m, ncol = p)
  RM <- RM[rowSums(abs(RM)) > 0, ]
  RM <- Matrix::Matrix(RM, sparse = TRUE)
  return(RM)
}

#' Gaussian Random Projection Matrix
#'
#' @param m goal dimension
#' @param p dimension of predictors to be projected
#'
#' @return an m x p matrix
#' @export
#'
#' @examples
#' generate_gaussian_rp(10, 20)
generate_gaussian_rp <- function(m, p) {
  vals <- rnorm(m * p)
  RM <- matrix(vals, nrow = m, ncol = p)
  RM <- Matrix::Matrix(RM, sparse = TRUE)
  return(RM)
}
