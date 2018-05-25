## test kenrl minque
library(MASS)
## library(microbenchmark)
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")

test <- function(N=100, P=200, r=100)
{
    X <- matrix(rnorm(N * P), N, P)

    e1 <- diag(N)
    p1 <- ply(X, coef0=1, degree=1)
    p2 <- ply(X, coef0=1, degree=2)
    p3 <- ply(X, coef0=1, degree=3)
    ga <- gau(X)

    ## allowed kernels
    V <- list(e1=e1, p1=p1)
    k <- length(V)

    ## contrast
    P <- rbind(diag(k), rep(1, k))

    ## true covariance
    S <- .1 * e1 + 1 * p1 + 2.0 * p2 + 1.5 * p3
    y <- mvrnorm(1, rep(0, N), S)

    r1 <- pkn.mnq.R(y, V, P, const=1, order=2)
    r1
}

perm <- function(n=3, order=3)
{
    r <- list(0:n)
    i <- 1
    while(i <= order)
    {
        . <- list(rep(tail(r, 1)[[1]], each=n))
        r <- c(r, .)
        i <- i + 1
    }
    do.call(cbind, r)
}
