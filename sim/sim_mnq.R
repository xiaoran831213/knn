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
    V <- knl[c('e', 'l', 'g')]

    ## true covariance
    X <- matrix(0, N, 1)
    S <- with(knl, .1 * e + 1. * l + 1. * g)
    y <- mvrnorm(1, rep(0, N), S)

    a0 <- list(
        a1=LMM_MINQUE_Solver(V, X, c(1, 0, 0))$A,
        a2=LMM_MINQUE_Solver(V, X, c(0, 1, 0))$A,
        a3=LMM_MINQUE_Solver(V, X, c(0, 0, 1))$A)
    a0 <- lapply(a0, round, 4)

    ## minque
    P <- diag(length(V))
    r1 <- knl.mnq.R(y, V, P, psd=FALSE)
    a1 <- lapply(r1$A, round, 4)
    
    r2 <- .Call('_knn_knl_mnq', PACKAGE = 'knn', as.matrix(y), V, P, 0)
    a2 <- lapply(r2$A, round, 4)
    
    for(i in seq_along(a0))
    {
        print(all.equal(a0[[i]], a1[[i]]))
        print(all.equal(a0[[i]], a2[[i]]))
    }

    list(a0=a0, r1=r1, r2=r2)
}

ts2 <- function(N=1000, P=2000, r=100)
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
    X <- matrix(0, N, 1)
    S <- with(knl, .1 * e + 1. * l + 1. * g)
    y <- mvrnorm(1, rep(0, N), S)

    a0 <- list(
        a1=Modified_MINQUE_Solver(V, X, c(1, 0, 0)),
        a2=Modified_MINQUE_Solver(V, X, c(0, 1, 0)),
        a3=Modified_MINQUE_Solver(V, X, c(0, 0, 1)))
    b0 <- unname(sapply(a0, function(a) t(y) %*% a %*% y))
    a0 <- lapply(a0, round, 5)

    ## minque
    P <- diag(length(V))
    r1 <- knl.mnq.R(y, V, P, psd=TRUE)
    a1 <- lapply(r1$A, round, 5)
    b1 <- unname(drop(r1$f))

    r2 <- .Call('_knn_knl_mnq', PACKAGE = 'knn', as.matrix(y), V, P, 1)
    a2 <- lapply(r2$A, round, 5)
    b2 <- unname(drop(r2$f))

    print(all.equal(b1, b2))
    list(a0=a0, r1=r1, r2=r2)
}

ts3 <- function(N=1000, P=2000, r=100)
{
    X <- matrix(rnorm(N * P), N, P)
    knl <- list(
        e=diag(N),
        l=ply(X, degree=1),
        q=ply(X, degree=2))
    ## g=gau(X))

    ## allowed kernels
    V <- knl[c('l')]

    ## true covariance
    X <- matrix(0, N, 1)
    S <- with(knl, .1 * e + 1. * l + 1. * q)
    y <- mvrnorm(1, rep(0, N), S)

    r1 <- knl.mnq(y, V, order=1, cpp=TRUE, psd=TRUE)
    r2 <- knl.mnq(y, V, order=1, cpp=FALSE, psd=TRUE)

    list(r1=r1, r2=r2)
}
