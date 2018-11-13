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
source("sim/mdl.R")
source("sim/mtd.R")
source("sim/gc1.R")
source("sim/grm.R")

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

ts2 <- function(N=200, P=500, r=10)
{
    xmx <- matrix(rnorm(N * P, 1), N, P)
    knl <- krn(xmx, ~LN1)               # working kernel
    fnl <- krn(xmx, ~GS1)               # function kernel
    
    ## true covariance
    W <- c(.5, 1, 2, 1, 1)[seq(1 + length(fnl))]
    C <- cmb(c(EPS(xmx), fnl), W)[[1]]  # function covariance
    y <- mvrnorm(1, rep(0, N), C)

    print('PDS=1, MINQUE')
    r1 <- knl.mnq(y, knl, NULL, itr=50, cpp=FALSE, psd=1)

    print('PDS=0, MINQUE')
    r0 <- knl.mnq(y, knl, NULL, itr=50, cpp=FALSE, psd=0)

    print('PDS=0, MINQUE')
    r2 <- rop.vcm(y, knl, NULL, cpp=FALSE)

    list(mn1=r1, mn0=r0, mle=r2, ref=W)
}

ts3 <- function(N=1000, P=2000, r=10)
{
    Z <- matrix(rpois(N * P, 2), N, P)

    ## fix effect to constitue the mean
    ## covariants
    X <- cbind(x1=rbinom(N, 1, .3), x2=rbinom(N, 1, .5), x3=rnorm(N))
    bts <- c(X1=1, X2=2, X3=1.5)
    M <- X %*% bts                      # mean
    ## X <- NULL
    ## M <- rep(0, N)
    
    ## random effect to constitute the covariance
    knl <- krn(Z, ~LN1)
    vcs <- c(LN1=2.0)
    eps <- c(EPS=0.5)
    S <- cmb(knl, vcs, drop=TRUE)       # Sigma

    ## response
    rsp <- mvrnorm(1, M, S) + rnorm(N, 0, sqrt(eps))

    ## working kernel
    use <- krn(Z, ~ID + LN1)

    ## mbk <- microbenchmark(
    ##     r1 <- .mn0(rsp, use, X=X)$vcs,
    ##     r2 <- .mnq(rsp, use, X=X)$vcs,
    ##     times=r)
    ## print(mbk)

    ## only once.
    r1 <- .mn0(rsp, use, X=X)$vcs
    r2 <- .mnq(rsp, use, X=X)$vcs
    r3 <- gct.rml(rsp, knl, qcvr=cbind(1, X))$par

    list(mnq.old=r1, mnq.new=r2, gct.rml=r3, ref=c(eps, vcs))
}
