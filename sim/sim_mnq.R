## test kenrl minque
library(MASS)
library(microbenchmark)
library(devtools)
devtools::load_all()
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("R/utl.R")
source("R/vcm.R")
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

ts3 <- function(N=1000, P=2000, bmk=1)
{
    Z.dvp <- matrix(rpois(N * P, 2), N, P)
    Z.evl <- matrix(rpois(N * P, 2), N, P)

    ## ------------------------------ generation ------------------------------ ##
    ## fix effect to constitue the mean
    ## intercept and covariants
    X.dvp <- cbind(X00=1, X01=rbinom(N, 1, .3), X02=rbinom(N, 1, .5), X03=rnorm(N))
    X.evl <- cbind(X00=1, X01=rbinom(N, 1, .3), X02=rbinom(N, 1, .5), X03=rnorm(N))
    bts <- c(X00=-0.6, X01=0.5, X02=1.0, X03=-1.0) # beta
    M.dvp <- X.dvp %*% bts                         # mean
    M.evl <- X.evl %*% bts                         # mean
    
    ## random effect to constitute the covariance
    K.dvp <- krn(Z.dvp, ~LN1)
    K.evl <- krn(Z.evl, ~LN1)
    eps <- c(EPS=1.5)
    vcs <- c(LN1=2.0)
    S.dvp <- cmb(K.dvp, vcs, drop=TRUE)       # Sigma
    S.evl <- cmb(K.evl, vcs, drop=TRUE)       # Sigma

    ## response
    y.dvp <- mvrnorm(1, M.dvp, S.dvp) + rnorm(N, 0, sqrt(eps))
    y.evl <- mvrnorm(1, M.evl, S.evl) + rnorm(N, 0, sqrt(eps))

    ## ------------------------------ working fit ------------------------------ ##
    ## working kernel
    X.use <- X.dvp # [, -1]
    K.use <- krn(Z.dvp, ~EPS + LN2)

    ## estimation of vcs
    if(bmk > 1)
    {
        bmk <- microbenchmark(
            vc1 <- .mn0(y.dvp, K.use, X=X.use)$vcs,
            vc2 <- .mnq(y.dvp, K.use, X=X.use)$vcs,
            times=bmk)
        print(bmk)
    }
    else
    {
        vc1 <- .mn0(y.dvp, K.use, X=X.use)$vcs
        vc2 <- .mnq(y.dvp, K.use, X=X.use)$vcs
    }
    md3 <- knl.mnq(y.dvp, K.use[-1], X.use, cpp=FALSE)
    md4 <- gct.rml(y.dvp, K.use[-1], X.use[, -1])
    
    ## estimation of covariance
    .v1 <- cmb(K.use, vc1, drop=TRUE)
    .v2 <- cmb(K.use, vc2, drop=TRUE)

    ## estimation of beta
    bt1 <- .gls(.v1, X.use, y.dvp)
    bt2 <- .gls(.v2, X.use, y.dvp)

    ## estimated parameters
    es1 <- list(par=c(bt1, vc1))
    es2 <- list(par=c(bt2, vc2))
    dvp <- list(es1=es1, es2=es2, es3=md3, es4=md4)
    par <- pars(dvp, ref=c(bts, eps, vcs))

    ## ------------------------------ testing ------------------------------ ##
    X.use <- X.evl # [, -1]
    K.use <- krn(Z.evl, ~EPS + LN2)
    pd1 <- vpd(y.evl, K=K.use[-1], es1$par, X.use)
    pd2 <- vpd(y.evl, K=K.use[-1], es2$par, X.use)
    pd3 <- vpd(y.evl, K=K.use[-1], md3$par, X.use)
    pd4 <- vpd(y.evl, K=K.use[-1], md4$par, X.use)
    evl <- list(ev1=pd1, ev2=pd2, ev3=pd3, ev4=pd4)
    
    ret <- list(par=par, dvp=dvp, evl=evl)
    print(ret$par)
    ret
}
