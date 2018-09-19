## test kenrl minque
library(MASS)
library(microbenchmark)
library(devtools)
devtools::load_all()
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("R/utl.R")
source("sim/utl.R")
source("sim/mn1.R")

## test for cpp minque versus R minque
ts1 <- function(N=2000, P=4000, r=5)
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
                         r2 <- .Call('knl_mnq', as.matrix(y), V, FALSE),
                         r3 <- .Call('egn_mnq', as.matrix(y), V), times=r)
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

ts2 <- function(N=100, P=50, frq=.1, r=10)
{
    X <- matrix(rnorm(N * P, 1), N, P)
    ## . <- sample(c(rep(TRUE, P*frq), rep(FALSE, P - P*frq)))
    ## F <- X[, .]

    knl <- krn(X, c(id, ki, s1))
    fnl <- krn(X, c(id, ki, s1))
    
    ## true covariance
    W <- c(eps=.5, vc1=1, vc2=2, vc3=1, vc4=2)[1:length(fnl)]
    y <- mvrnorm(1, rep(0, N), cmb(fnl, W)[[1]])

    r1 <- knl.mnq.R(y, knl, X=NULL)$vcs
    r2 <- itr.mnq(y, knl, X=NULL)$par
    r3 <- round(drop(rop.vcm(y, knl[-1])$par), 4)

    list(mnq.old=r1, mnq.new=r2, mle=r3, ref=W)
}

ts3 <- function(N=100, P=50, frq=.1, r=10)
{
    X1 <- matrix(rpois(N * P, 2), N, P)
    X2 <- matrix(rnorm(N * P, 2), N, P)
    K1 <- krn(X1, p1)[[1]]
    K2 <- krn(X2, ga)[[1]]
    id <- diag(N)

    knl <- list(id, K1)
    use <- list(id, K1, K2)
    
    ## true covariance
    vcs <- c(eps=1, vc1=1, vc2=1)[1:length(knl)]
    cmx <- cmb(knl, vcs)[[1]]
    y <- mvrnorm(1, rep(0, N), cmx)

    r1 <- knl.mnq.R(y, use, X=NULL)$vcs
    r2 <- itr.mnq(y, use, X=NULL)$par
    r3 <- round(drop(rop.vcm(y, use[-1])$par), 4)

    list(mnq.old=r1, mnq.new=r2, mle=r3, ref=vcs)
}
