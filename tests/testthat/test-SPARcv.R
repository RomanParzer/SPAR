
# similar to spar

test_that("Results has right class", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y)
  expect_equal(class(spar_res),"spar.cv")
})

test_that("Coef returns vector of correct length", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x,y)
  sparcoef <- coef(spar_res)
  expect_equal(length(sparcoef$beta),30)
})

test_that("Coef is more sparse for 1se rule", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x,y)
  sparcoef <- coef(spar_res,opt_par = "best")
  sparcoef2 <- coef(spar_res,opt_par = "1se")
  expect_equal(all(which(sparcoef$beta==0) %in% which(sparcoef2$beta==0)),TRUE)
})

test_that("Validated lambda values are same as the ones for initial SPAR fit", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y,nummods=c(10,15))
  expect_equal(unique(spar_res$val_sum$lam),as.numeric(spar_res$lambdas))
})

test_that("Validated nummod values are same as the ones for initial SPAR fit", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y,nummods=c(10,15))
  expect_equal(unique(spar_res$val_sum$nummod),as.numeric(spar_res$nummods))
})

# Tests expecting errors

test_that("Get errors for input x not data.frame or matrix", {
  x <- list("1"=1:10,"2"=(-1)^(1:12),"3"=rnorm(12),"4"=rnorm(12),"5"=runif(12),"6"=runif(12))
  y <- rnorm(6)
  expect_error(spar.cv(x,y,nfolds=2))
})

test_that("Get errors for input x non-numeric data.frame", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  x$X1[1] <- "a"
  y <- rnorm(10)
  expect_error(spar.cv(x,y))
})

test_that("Get errors for categorical input y", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- factor(rep(c("a","b"),5))
  expect_error(spar.cv(x,y))
})

test_that("Get errors for prediction when xnew has wrong dimensions", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar.cv(x,y)
  xnew <- example_data$xtest
  expect_error(predict(spar_res,xnew=xnew[,-1]))
})

