
test_that("Results has right class", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar(x,y)
  expect_equal(class(spar_res),"spar")
})

test_that("Coef returns vector of correct length", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar(x,y)
  sparcoef <- coef(spar_res)
  expect_equal(length(sparcoef$beta),30)
})

test_that("Coef is more sparse for higher threshold", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar(x,y)
  sparcoef <- coef(spar_res)
  sparcoef2 <- coef(spar_res,lambda = spar_res$val_res$lam[which(spar_res$val_res$lam==sparcoef$lambda)+1])
  expect_equal(all(which(sparcoef$beta==0) %in% which(sparcoef2$beta==0)),TRUE)
})

test_that("Returned coef and preds are correct for fixed screening and projections", {
  x <- example_data$x
  y <- example_data$y
  xnew <- example_data$xtest

  m <- 2*floor(log(ncol(x)))
  nsc <- 2*nrow(x)
  RP1 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP1@i <- as.integer(c(rep(1:m,each=nsc%/%m),rep(m,nsc%%m))-1)
  RP1@p <- 0:nsc
  RP1@x <- (-1)^(1:nsc)

  m <- floor(nrow(x)/2)
  RP2 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP2@i <- as.integer(c(rep(m:1,each=nsc%/%m),rep(1,nsc%%m))-1)
  RP2@p <- 0:nsc
  RP2@x <- (-1)^(1:nsc)

  spar_res <- spar(x,y,nummods=c(2),inds = list(1:(2*nrow(x)),500+1:(2*nrow(x))),RPMs = list(RP1,RP2))
  sparcoef <- coef(spar_res)
  pred <- predict(spar_res,xnew=xnew)

  expect_equal(sparcoef$lambda,0.002324331,tolerance = 1e-6)
  expect_equal(sparcoef$beta[53],0)
  expect_equal(sparcoef$beta[1],0.1384109,tolerance = 1e-6)
  expect_equal(pred[1],20.78687,tolerance = 1e-5)
})

test_that("Returned coef and preds are correct for fixed screening and projections for binomial(logit)", {
  x <- example_data$x
  y <- round(1/(1+exp(-example_data$y)))
  xnew <- example_data$xtest

  m <- 2*floor(log(ncol(x)))
  nsc <- 2*nrow(x)
  RP1 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP1@i <- as.integer(c(rep(1:m,each=nsc%/%m),rep(m,nsc%%m))-1)
  RP1@p <- 0:nsc
  RP1@x <- (-1)^(1:nsc)

  m <- floor(nrow(x)/2)
  RP2 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP2@i <- as.integer(c(rep(m:1,each=nsc%/%m),rep(1,nsc%%m))-1)
  RP2@p <- 0:nsc
  RP2@x <- (-1)^(1:nsc)

  spar_res <- spar(x,y,nummods=c(2),inds = list(1:(2*nrow(x)),500+1:(2*nrow(x))),RPMs = list(RP1,RP2),family = binomial(logit))
  sparcoef <- coef(spar_res)
  pred <- predict(spar_res,xnew=xnew)

  expect_equal(sparcoef$lambda,0.009781171,tolerance = 1e-6)
  expect_equal(sparcoef$beta[11],0)
  expect_equal(sparcoef$beta[1],0.04584074,tolerance = 1e-6)
  expect_equal(pred[1],0.9675592,tolerance = 1e-5)
})

test_that("Columns with zero sd get ceofficient 0", {
  x <- example_data$x
  x[,c(1,11,111)] <- 2
  y <- example_data$y

  spar_res <- spar(x,y)
  sparcoef <- coef(spar_res)
  expect_equal(sparcoef$beta[c(1,11,111)],c(0,0,0))
})

# Tests expecting errors

test_that("Get errors for data.frame input x", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  expect_error(spar(x,y))
})

test_that("Get errors for categorical input y", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- factor(rep(c("a","b"),5))
  expect_error(spar(x,y))
})

test_that("Get errors for mslow > msup", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  expect_error(spar(x,y,mslow = 15,msup = 10))
})

test_that("Get errors for msup > nscreen", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  expect_error(spar(x,y,nscreen=18,msup = 20))
})

test_that("Get errors for to small length of inds and RPMs lists", {
  x <- example_data$x
  y <- example_data$y

  m <- 2*floor(log(ncol(x)))
  nsc <- 2*nrow(x)
  RP1 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP1@i <- as.integer(c(rep(1:m,each=nsc%/%m),rep(m,nsc%%m))-1)
  RP1@p <- 0:nsc
  RP1@x <- (-1)^(1:nsc)

  m <- floor(nrow(x)/2)
  RP2 <- Matrix::Matrix(c(0),nrow=m,ncol=nsc,sparse=TRUE)
  RP2@i <- as.integer(c(rep(m:1,each=nsc%/%m),rep(1,nsc%%m))-1)
  RP2@p <- 0:nsc
  RP2@x <- (-1)^(1:nsc)

  expect_error(spar(x,y,nummods=c(3),inds = list(1:(2*nrow(x)),500+1:(2*nrow(x))),RPMs = list(RP1,RP2)))
})

test_that("Get errors for prediction when xnew has wrong dimensions", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar(x,y)
  xnew <- example_data$xtest
  expect_error(predict(spar_res,xnew=xnew[,-1]))
})

test_that("Get errors for classification validation measure for non-binomial family", {
  x <- example_data$x
  y <- example_data$y
  expect_error(spar(x,y,type.measure = "1-auc"))
})


