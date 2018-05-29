## test kenrl minque
library(MASS)
## library(microbenchmark)
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")

ts1 <- function(N=1000, P=2000, r=100)
{
    X <- matrix(rnorm(N * P), N, P)
    knl <- list(
        e=diag(N),
        l=ply(X, degree=1),
        q=ply(X, degree=2),
        g=gau(X))

    ## allowed kernels
    V <- knl[c('e', 'l', 'g')]

    ## true covariance
    S <- with(knl, .1 * e + 1. * l, 1. * g)
    y <- mvrnorm(1, rep(0, N), S)

    ## minque
    r1 <- knl.mnq.R(y, V)
    r1
}

ts2 <- function(N=1000, P=2000, r=100)
{
    ## sumulated X, y, and covariance
    X <- matrix(rnorm(N * P), N, P)
    knl <- list(
        e=diag(N),
        l=ply(X, degree=1),
        q=ply(X, degree=2),
        g=gau(X, sigma=1.))
    
    S <- with(knl, .1 * e + 1. * g)
    y <- mvrnorm(1, rep(0, N), S)
    
    ## allowed kernels
    V <- knl[c('l')]

    r1 <- knl.mnq(y, V, cpp=TRUE, order=1)
    r2 <- knl.mnq(y, V, cpp=TRUE, order=2)
    r3 <- knl.mnq(y, V, cpp=TRUE, order=3)
    
    list(r1=r1, r2=r2, r3=r3)
}
