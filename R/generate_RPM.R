
#' Generate CW Random Projection Matrix (Definition 2.2. in Parzer, Vana_Gür, Filzmoser 2023)
#'
#' @param m goal dimension
#' @param p initial dimension
#' @param coef p-vector of diagonal elements, random +/-1 by default
#' @returns sparse mxp matrix of class "dgCMatrix" (Matrix package)
#' @keywords internal
generate_RPM <- function(m,
                         p,
                         coef=sample(c(-1,1),p,replace = TRUE)) {
  goal_dims <- sample(1:m,p,replace = TRUE)
  counter <- 0
  # remove zero rows
  for (goal_dim in 1:m) {
    if (sum(goal_dims==(goal_dim-counter))==0) {
      goal_dims[goal_dims>goal_dim-counter] <- goal_dims[goal_dims>goal_dim-counter]-1
      counter <- counter + 1
    }
  }

  RM <- Matrix::Matrix(c(0),nrow=m-counter,ncol=p,sparse=TRUE)
  RM@i <- as.integer(goal_dims-1)
  RM@p <- 0:p
  RM@x <- coef

  return(RM)
}
