## test kenrl minque
library(MASS)
## library(microbenchmark)
library(devtools)
devtools::load_all()
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")

## test for cpp minque versus R minque
ts1 <- function(N=1000, P=2000, r=100)
{
    X <- matrix(rnorm(N * P), N, P)
    knl <- list(
        e=diag(N),
        l=ply(X, degree=1),
        q=ply(X, degree=2),
        g=gau(X))

    ## allowed kernels
    V <- knl[c('e', 'l')]

    ## true covariance
    X <- matrix(0, N, 1)
    S <- with(knl, .1 * e + 1. * g)
    y <- mvrnorm(1, rep(0, N), S)

    ## minque
    r1 <- knl.mnq.R(y, V, X=NULL)
    a1 <- lapply(r1$A, round, 4)

    r2 <- .Call('knl_mnq', as.matrix(y), V)
    a2 <- lapply(r2$A, round, 4)

    for(i in seq_along(a1))
        print(all.equal(a1[[i]], a2[[i]]))
    
    list(r1=r1, r2=r2)
}
