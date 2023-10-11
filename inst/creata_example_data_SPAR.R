# # # create illustration data set from Example 2.1 in Parzer, Vana-Guer, Filzmoser (2023)
n <- 200
p <- 2000
ntest <- 100
a <- 100
snr <- 10
rho <- 0.5
mu <- 1

beta <- numeric(p)
set.seed(1234)
beta[1:a] <- sample(c(-3:3)[-4],a,replace = TRUE)
x <- sqrt(rho)*matrix(rep(rnorm((n+ntest),0,1),p),n+ntest,p) + sqrt(1-rho)*matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
bSb <- rho*sum(beta)^2 + (1-rho)*sum(beta^2)
sigma2 <- bSb/snr
y <- mu + x%*%beta + rnorm(n+ntest,0,sqrt(sigma2))

xtest <- x[-(1:n),]
x <- x[1:n,]
ytest <- y[-(1:n)]
y <- y[1:n]
example_data <- list(x=x,y=y,xtest=xtest,ytest=ytest,mu=mu,beta=beta,sigma2=sigma2)
save(example_data,file="../data/example_data.rda")
