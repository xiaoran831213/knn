## test kenrl minque
library(MASS)
library(microbenchmark)
library(devtools)
devtools::load_all()
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("sim/utl.R")

## test for cpp minque versus R minque
ts1 <- function(N=1000, P=2000, r=10)
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
    mb <- microbenchmark(r1 <- knl.mnq.R(y, V, X=NULL),
                         r2 <- .Call('knl_mnq', as.matrix(y), V), times=r)
    print(mb)

    ## for(i in seq_along(length(V)))
    ## {
    ##     print(all.equal(r1$A[[i]], r2$A[[i]]))
    ## }
    ## print(all.equal(r1$C, r2$C))
    ## print(all.equal(r1$se, drop(r2$se)))
    ## print(all.equal(r1$s2, drop(r2$s2)))

    ## list(r1=r1, r2=r2)
}

ts2 <- function(N=6, P=3, r=10)
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
    S <- with(knl, .01 * e + 1. * g)
    y <- mvrnorm(1, rep(0, N), S)

    mb <- microbenchmark(
        ## a <- knl.mnq.R(y, V, X=NULL),
        b <- .Call('knl_mnq', as.matrix(y), V),
        c <- .Call('egn_mnq', as.matrix(y), V))
    print(mb)

    list(
        #3 rln=a,
        am1=b, am2=c)
}
