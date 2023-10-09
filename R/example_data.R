#' Example dataset to illustrate spar functionalities
#'
#' High-dimensional regression dataset from Example 2.1 in Parzer, Vana-Guer, Filzmoser (2023) with n=200, p=2000, ntest=100,
#' with true parameters such that y = mu + x%*%beta + rnorm(n,0,sqrt(sigma2))
#'
#' @format A list with objects:
#' \describe{
#'   \item{x}{n x p predictor matrix from multivariate gaussian distribution}
#'   \item{y}{n-vector of responses}
#'   \item{xtest}{ntest x p predictor matrix independent of x}
#'   \item{ytest}{ntest-vector of test response observations}
#'   \item{mu}{true scalar mean used to generate responses}
#'   \item{beta}{true p-vector of coefficients}
#'   \item{sigma2}{true noise variance}
#' }
"example_data"
