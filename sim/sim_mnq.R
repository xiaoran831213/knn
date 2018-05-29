## test kenrl minque
library(MASS)
## library(microbenchmark)
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")

test <- function(N=100, P=200, r=100)
{
    X <- matrix(rnorm(N * P), N, P)

    knl <- list(
        e=diag(N),
        l=ply(X, coef0=1, degree=1),
        q=ply(X, coef0=1, degree=2),
        c=ply(X, coef0=1, degree=3),
        g=gau(X))

    ## allowed kernels
    V <- knl[c('l')]
    k <- length(V) + 1

    ## contrast
    P <- rbind(diag(k), rep(1, k))

    ## true covariance
    S <- with(knl, .1 * e + 1 * l + 2.0 * q + 1.5 * c)
    y <- mvrnorm(1, rep(0, N), S)

    r1 <- pkn.mnq.R(y, V, const=1, order=2)
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
